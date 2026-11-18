// Adapted from deep_autolabel memory_mode/spatial_sliding/window.h.
#pragma once
#include <Eigen/Core>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>

namespace flivo {

// World-aligned cube shared by LIO and VIO, independent of their voxel sizes.
struct SpatialWindowBounds {
    Eigen::Vector3d center;
    double half_extent_m;

    // Retain boundary voxels whole, including their tree/plane statistics.
    // Do not narrow the existing int64_t keys to int or add to them in int64_t.
    bool IntersectsVoxel(int64_t x, int64_t y, int64_t z, double size) const {
        const int64_t key[3] = {x, y, z};
        for (int axis = 0; axis < 3; ++axis) {
            const long double lo = static_cast<long double>(key[axis]) * size;
            const long double hi = (static_cast<long double>(key[axis]) + 1) * size;
            if (hi < static_cast<long double>(center[axis]) - half_extent_m ||
                lo > static_cast<long double>(center[axis]) + half_extent_m)
                return false;
        }
        return true;
    }
};

class SpatialSlidingWindow {
 public:
    SpatialSlidingWindow(double half_extent_m, double trigger_m)
        : half_extent_m_(half_extent_m), trigger_m_(trigger_m) {
        if (!std::isfinite(half_extent_m) || !std::isfinite(trigger_m) ||
            half_extent_m <= 0 || trigger_m <= 0 || trigger_m >= half_extent_m)
            throw std::invalid_argument("invalid spatial sliding window dimensions");
    }

    bool Update(const Eigen::Vector3d& position) {
        if (!position.allFinite())
            throw std::invalid_argument("nonfinite spatial sliding window position");
        if (bounds_ && (position - bounds_->center).norm() < trigger_m_)
            return false;
        bounds_ = SpatialWindowBounds{position, half_extent_m_};
        return true;
    }

    void Reset() { bounds_.reset(); }
    const std::optional<SpatialWindowBounds>& bounds() const { return bounds_; }

 private:
    double half_extent_m_, trigger_m_;
    std::optional<SpatialWindowBounds> bounds_;
};

}  // namespace flivo
