#pragma once
// Adapted from deep_autolabel frontend/fast-livo2/page/image_backing.h.
// Native FAST-LIVO2 references one grayscale image, rather than a pyramid.
#include <opencv2/core.hpp>
#include <cerrno>
#include <filesystem>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>
#include <sys/mman.h>
#include <unistd.h>

namespace flivo {
struct MappedImage {
    cv::Mat pixels;
    std::shared_ptr<void> owner;
    size_t bytes = 0;
};

// Reference observations own the mapping, even after the store is destroyed.
// MADV_DONTNEED releases PTEs; later reads fault in the same file-backed bytes.
class ImageBackingStore {
    struct Mapping {
        void* data = MAP_FAILED;
        size_t bytes = 0;
        ~Mapping() { if (data != MAP_FAILED) munmap(data, bytes); }
    };
public:
    struct Stats {
        size_t images = 0, mapped_bytes = 0, working_bytes = 0;
        uint64_t drop_calls = 0;
    } stats;

    explicit ImageBackingStore(const std::filesystem::path& root) : root_(root) {
        std::filesystem::create_directories(root_);
    }
    std::shared_ptr<const MappedImage> Store(const cv::Mat& source) {
        if (source.empty()) throw std::runtime_error("cannot back an empty reference image");
        std::string name = (root_ / "image.XXXXXX").string();
        std::vector<char> path(name.begin(), name.end()); path.push_back(0);
        const int fd = mkstemp(path.data());
        if (fd < 0) throw std::runtime_error("cannot create reference image backing");
        if (unlink(path.data()) != 0) { close(fd); throw std::runtime_error("cannot unlink reference image backing"); }
        auto backing = std::make_shared<Mapping>();
        try {
            for (int row = 0; row < source.rows; ++row) {
                const char* data = reinterpret_cast<const char*>(source.ptr(row));
                size_t remaining = source.cols * source.elemSize(); backing->bytes += remaining;
                while (remaining) {
                    const auto n = write(fd, data, remaining);
                    if (n < 0 && errno == EINTR) continue;
                    if (n <= 0) throw std::runtime_error("reference image backing write failed");
                    data += n; remaining -= n;
                }
            }
            backing->data = mmap(nullptr, backing->bytes, PROT_READ, MAP_SHARED, fd, 0);
            if (backing->data == MAP_FAILED) throw std::runtime_error("reference image mmap failed");
        } catch (...) { close(fd); throw; }
        close(fd);
        auto out = std::make_shared<MappedImage>();
        out->pixels = cv::Mat(source.rows, source.cols, source.type(), backing->data);
        out->owner = backing; out->bytes = backing->bytes;
        records_.push_back(out);
        return out;
    }
    void Trim(const std::set<const MappedImage*>& working) {
        stats.images = stats.mapped_bytes = stats.working_bytes = 0;
        for (auto it = records_.begin(); it != records_.end();) {
            const auto image = it->lock();
            if (!image) { it = records_.erase(it); continue; }
            ++stats.images; stats.mapped_bytes += image->bytes;
            if (working.count(image.get())) stats.working_bytes += image->bytes;
            else {
                if (madvise(image->pixels.data, image->bytes, MADV_DONTNEED) != 0)
                    throw std::runtime_error("reference image page release failed");
                ++stats.drop_calls;
            }
            ++it;
        }
    }
private:
    std::filesystem::path root_;
    std::vector<std::weak_ptr<const MappedImage>> records_;
};
} // namespace flivo
