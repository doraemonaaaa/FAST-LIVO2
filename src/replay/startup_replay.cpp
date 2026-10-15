/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include "replay/startup_replay.h"
#include "LIVMapper.h"
#include <iomanip>
#include <map>
#include <stdexcept>

void StartupReplay::readParameters(ros::NodeHandle &nh, bool imu_en)
{
  nh.param<bool>("imu/initialization_replay", enabled_, false);
  int initialization_replay_max_mb;
  nh.param<int>("imu/initialization_replay_max_mb", initialization_replay_max_mb, 2048);
  if (enabled_ && (!imu_en || initialization_replay_max_mb <= 0))
    throw std::runtime_error("Initialization replay requires enabled IMU and a positive cache limit");
  limit_ = size_t(initialization_replay_max_mb) * 1024 * 1024;
}

void StartupReplay::initializeLog(const std::string &path)
{
  if (enabled_)
  {
    log_.open(path);
    if (!log_) throw std::runtime_error("Cannot create initialization replay log");
    log_ << "event,begin,end,groups,bytes\n" << std::flush;
  }
}

void StartupReplay::cacheImu(const sensor_msgs::Imu::ConstPtr &imu)
{
  if (!pending()) return;
  if (bytes_ > limit_ - sizeof(sensor_msgs::Imu))
    throw std::runtime_error("Initialization replay IMU cache limit exceeded");
  imu_samples_.push_back(imu);
  bytes_ += sizeof(sensor_msgs::Imu);
}

void StartupReplay::cacheMeasurement(const LidarMeasureGroup &group)
{
  if (!pending()) return;
  size_t bytes = sizeof(LidarMeasureGroup);
  for (const auto &cloud : {group.lidar, group.pcl_proc_cur, group.pcl_proc_next})
    if (cloud) bytes += cloud->size() * sizeof(PointType);
  for (const auto &m : group.measures)
  {
    bytes += m.imu.size() * sizeof(sensor_msgs::Imu);
    for (const auto &image : m.images) bytes += image.pixels.total() * image.pixels.elemSize();
  }
  if (bytes > limit_ || bytes_ > limit_ - bytes)
    throw std::runtime_error("Initialization replay cache limit exceeded; no startup data was silently discarded");
  LidarMeasureGroup saved = group;
  if (group.lidar) saved.lidar.reset(new PointCloudXYZI(*group.lidar));
  if (group.pcl_proc_cur) saved.pcl_proc_cur.reset(new PointCloudXYZI(*group.pcl_proc_cur));
  if (group.pcl_proc_next) saved.pcl_proc_next.reset(new PointCloudXYZI(*group.pcl_proc_next));
  for (auto &m : saved.measures) m.cloneImages();
  measurements_.push_back(std::move(saved));
  bytes_ += bytes;
}

void StartupReplay::replay(LIVMapper &mapper)
{
  if (measurements_.empty()) throw std::runtime_error("Empty initialization replay cache");
  // Keep synchronization's live pending scan/image state separate from replay.
  LidarMeasureGroup live = mapper.LidarMeasures;
  const double begin = measurements_.front().last_lio_update_time;
  const double end = live.last_lio_update_time;
  std::map<uint64_t, sensor_msgs::Imu::ConstPtr> ordered;
  // Include IMU received before the first LiDAR and samples omitted by slicing.
  for (const auto &imu : imu_samples_)
    if (imu->header.stamp.toSec() <= end) ordered.emplace(imu->header.stamp.toNSec(), imu);
  for (const auto &group : measurements_)
    for (const auto &m : group.measures)
      for (const auto &imu : m.imu) ordered.emplace(imu->header.stamp.toNSec(), imu);
  deque<sensor_msgs::Imu::ConstPtr> samples;
  for (const auto &item : ordered)
  {
    if (item.second->header.stamp.toSec() <= begin) samples.clear();
    samples.push_back(item.second);
  }
  mapper.p_imu->rewindForStartupReplay(mapper._state, end, begin, samples);
  mapper.vio_manager->clearStartupMap();
  for (auto &item : mapper.voxelmap_manager->voxel_map_) delete item.second;
  mapper.voxelmap_manager->voxel_map_.clear(); mapper.voxelmap_manager->ptpl_list_.clear();
  mapper.voxelmap_manager->pv_list_.clear(); mapper._pv_list.clear();
  mapper.voxelmap_manager->last_slide_position = mapper._state.pos_end;
  mapper.lidar_map_inited = false;
  for (auto &cloud : {mapper.feats_undistort, mapper.feats_down_body, mapper.feats_down_world,
      mapper.pcl_w_wait_pub, mapper.pcl_b_wait_pub, mapper.pcl_wait_pub, mapper.pcl_b_wait_color}) cloud->clear();
  mapper.pcl_wait_save->clear(); mapper.pcl_wait_save_intensity->clear();
  mapper.path.poses.clear();
  replay_started_ = true;
  replaying_ = true;
  log_ << std::setprecision(17) << "begin," << begin << ',' << end << ','
      << measurements_.size() << ',' << bytes_ << '\n' << std::flush;
  ROS_WARN("[startup replay] replaying %zu groups from %.9f to %.9f; cache %.1f MiB",
      measurements_.size(), begin, end, bytes_ / 1048576.);
  size_t count = 0;
  while (!measurements_.empty())
  {
    mapper.LidarMeasures = std::move(measurements_.front());
    measurements_.pop_front();
    mapper.processImu();
    mapper.stateEstimationAndMapping();
    ++count;
    log_ << std::setprecision(17) << "group,"
        << mapper.LidarMeasures.lidar_frame_beg_time << ',' << mapper.LidarMeasures.last_lio_update_time
        << ',' << count << ",0\n";
  }
  mapper.LidarMeasures = live;
  replaying_ = false;
  log_ << std::setprecision(17) << "complete," << begin << ',' << end
      << ',' << count << ',' << bytes_ << '\n' << std::flush;
  bytes_ = 0;
  imu_samples_.clear();
}

void StartupReplay::finish()
{
  log_.flush();
  if (pending() && !measurements_.empty())
    throw std::runtime_error("Recording ended before IMU initialization; startup replay is incomplete");
}
