#include "auxiliary_lidar_export.h"
#include "IMU_Processing.h"
#include "camera_imu_motion.h"
#include "imu_support.h"
#include "preprocess.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <ros/subscriber.h>
#include <xmlrpcpp/XmlRpcValue.h>

namespace
{
constexpr double kMaximumPoseAge = 0.2;
constexpr double kPoseHistorySeconds = 1.0;
constexpr double kImuHistorySeconds = 5.0;
constexpr size_t kPendingScans = 32;

struct Scan
{
  uint64_t stamp;
  double begin, end;
  PointCloudXYZI::Ptr points;
};
}

struct AuxiliaryLidarExport::Stream
{
  Stream(ros::NodeHandle &nh, const std::string &name, const std::filesystem::path &root,
         double blind_range, bool save_scans) : directory(root / name)
  {
    const auto key = "auxiliary_lidar_export/" + name + "/";
    std::string topic;
    int type, scan_lines;
    std::vector<double> rotation, translation;
    if (!nh.getParam(key + "topic", topic) || !nh.getParam(key + "lidar_type", type) ||
        !nh.getParam(key + "scan_lines", scan_lines) || scan_lines < 1 || scan_lines > 128 ||
        !nh.getParam(key + "body_from_lidar/rotation", rotation) || rotation.size() != 9 ||
        !nh.getParam(key + "body_from_lidar/translation_m", translation) || translation.size() != 3 ||
        (type != ROBOSENSE && type != PANDAR128))
      throw std::runtime_error("Invalid auxiliary LiDAR profile: " + name);

    rotation_to_imu << MAT_FROM_ARRAY(rotation);
    translation_to_imu << VEC_FROM_ARRAY(translation);
    if (!rotation_to_imu.allFinite() || !translation_to_imu.allFinite())
      throw std::runtime_error("Nonfinite auxiliary LiDAR calibration: " + name);
    preprocess.set(false, type, blind_range, 1);
    preprocess.N_SCANS = scan_lines;
    preprocess.blind_sqr = blind_range * blind_range;
    if (save_scans)
    {
      std::filesystem::create_directories(directory);
      const auto path = directory / "index.csv";
      if (std::filesystem::exists(path))
        throw std::runtime_error("Auxiliary LiDAR output already exists: " + path.string());
      index.open(path);
      index.imbue(std::locale::classic());
      index << "scan_start_s,scan_end_s,points,reference_timestamp,status,reason,pcd\n" << std::flush;
      if (!index) throw std::runtime_error("Cannot create " + path.string());
    }
    subscriber = nh.subscribe(topic, kPendingScans * 2, &Stream::receive, this);
  }

  void record(const Scan &scan, const Reference *ref, const char *reason,
              const std::string &filename = "")
  {
    index << std::fixed << std::setprecision(9) << scan.begin << ',' << scan.end
          << ',' << (scan.points ? scan.points->size() : 0) << ',';
    if (ref) index << std::setprecision(6) << ref->time;
    index << ',' << (filename.empty() ? "dropped" : "accepted") << ',' << reason
          << ',' << filename << '\n' << std::flush;
    if (!index) throw std::runtime_error("Cannot write auxiliary LiDAR scan index");
  }

  void receive(const sensor_msgs::PointCloud2::ConstPtr &msg)
  {
    // During initialization preparation, rosbag waits for this subscriber
    // alongside the estimator inputs. The second pass rereads the scans.
    if (!index.is_open()) return;
    const double begin = msg->header.stamp.toSec();
    Scan scan{msg->header.stamp.toNSec(), begin, begin, {}};
    if (begin <= last_begin)
    {
      record(scan, nullptr, "nonmonotonic_scan");
      return;
    }
    last_begin = begin;
    if (!msg->width || !msg->height || msg->data.empty())
    {
      record(scan, nullptr, "empty_scan");
      return;
    }

    scan.points.reset(new PointCloudXYZI());
    preprocess.process(msg, scan.points);
    if (!scan.points->empty())
      scan.end += .001 * scan.points->points.back().curvature;
    if (scan.points->empty() || !std::isfinite(scan.end) || scan.end <= begin ||
        scan.points->points.front().curvature < 0)
    {
      record(scan, nullptr, "invalid_scan_interval");
      return;
    }
    if (scans.size() == kPendingScans)
    {
      record(scans.front(), nullptr, "buffer_overflow");
      scans.pop_front();
    }
    scans.push_back(std::move(scan));
  }

  std::filesystem::path directory;
  std::ofstream index;
  Preprocess preprocess;
  M3D rotation_to_imu;
  V3D translation_to_imu;
  double last_begin = -1;
  std::deque<Scan> scans;
  ros::Subscriber subscriber;
};

AuxiliaryLidarExport::AuxiliaryLidarExport(ros::NodeHandle &nh, const std::string &directory,
                                         double blind_range, bool save_scans)
{
  XmlRpc::XmlRpcValue sources;
  if (!nh.getParam("auxiliary_lidar_export", sources) || sources.getType() != XmlRpc::XmlRpcValue::TypeStruct)
    throw std::runtime_error("auxiliary_lidar_export must be a map of sensor profiles");
  for (const auto &source : sources)
    streams_.push_back(std::make_unique<Stream>(nh, source.first, directory, blind_range, save_scans));
}

AuxiliaryLidarExport::~AuxiliaryLidarExport() = default;

void AuxiliaryLidarExport::observeImu(const sensor_msgs::Imu::ConstPtr &sample)
{
  samples_.push_back(sample);
  const double cutoff = imu_support::Time(sample) - kImuHistorySeconds;
  while (samples_.size() > 2 && imu_support::Time(samples_[1]) <= cutoff)
    samples_.pop_front();
}

void AuxiliaryLidarExport::observePose(double time, const StatesGroup &state)
{
  if (!references_.empty() && time <= references_.back().time)
    throw std::runtime_error("Auxiliary LiDAR export requires increasing main pose timestamps");
  references_.push_back({time, state});
  while (references_.size() > 2 && references_[1].time < time - kPoseHistorySeconds)
    references_.pop_front();
}

void AuxiliaryLidarExport::drain(const ImuProcess &imu, double estimator_time, bool finished)
{
  for (auto &stream : streams_)
  {
    while (!stream->scans.empty())
    {
      const Scan &scan = stream->scans.front();
      auto next = std::lower_bound(references_.begin(), references_.end(), scan.end,
          [](const Reference &ref, double time) { return ref.time < time; });
      // Wait for the next possible pose until the estimator has progressed
      // past the pairing window. At EOF the latest pose can still serve.
      if (next == references_.end() && !finished && estimator_time < scan.end + kMaximumPoseAge)
        break;
      auto nearest = next;
      if (next != references_.begin())
      {
        auto previous = std::prev(next);
        if (next == references_.end() || scan.end - previous->time <= next->time - scan.end)
          nearest = previous;
      }
      const Reference *ref = nearest == references_.end() ? nullptr : &*nearest;
      if (!ref || std::abs(ref->time - scan.end) > kMaximumPoseAge)
      {
        stream->record(scan, nullptr, "no_nearby_main_pose");
      }
      else
      {
        const auto coverage = imu_support::Check(samples_, std::min(scan.begin, ref->time),
                                                 std::max(scan.end, ref->time));
        if (coverage == imu_support::Coverage::MissingEnd && !finished) break;
        if (coverage != imu_support::Coverage::Covered)
        {
          stream->record(scan, ref, imu_support::Reason(coverage));
        }
        else
        {
          const flivo::CameraImuMotion motion(samples_, ref->time, scan.end, G_m_s2 / imu.IMU_mean_acc_norm);
          const StatesGroup at_end = motion.evaluate(ref->state);
          PointCloudXYZI::Ptr deskewed(new PointCloudXYZI());
          imu.deskewWholeScan(*scan.points, scan.begin, scan.end, at_end, samples_,
                              stream->rotation_to_imu, stream->translation_to_imu, deskewed);
          if (deskewed->size() != scan.points->size())
            throw std::runtime_error("Auxiliary whole-scan deskew changed point count");

          const M3D rotation = ref->state.rot_end.transpose() * at_end.rot_end;
          const V3D translation = ref->state.rot_end.transpose() * (at_end.pos_end - ref->state.pos_end);
          const M3D rotation_to_reference = rotation * stream->rotation_to_imu;
          const V3D translation_to_reference = rotation * stream->translation_to_imu + translation;
          for (auto &point : deskewed->points)
          {
            const V3D xyz = rotation_to_reference * V3D(point.x, point.y, point.z) + translation_to_reference;
            point.getVector3fMap() = xyz.cast<float>();
          }
          deskewed->width = deskewed->size();
          deskewed->height = 1;
          const std::string filename = std::to_string(scan.stamp) + ".pcd";
          if (pcl::PCDWriter().writeBinary((stream->directory / filename).string(), *deskewed) != 0)
            throw std::runtime_error("Cannot save auxiliary LiDAR scan: " + filename);
          stream->record(scan, ref, "covered", filename);
        }
      }
      stream->scans.pop_front();
    }
  }
}
