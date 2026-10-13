/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#ifndef FAST_LIVO_STARTUP_REPLAY_H
#define FAST_LIVO_STARTUP_REPLAY_H

#include "common_lib.h"
#include <deque>
#include <fstream>
#include <ros/node_handle.h>

class LIVMapper;

// Owns startup data until initialization and map admission allow chronological
// replay. Sensor synchronization and normal LIO/VIO updates remain in LIVMapper.
class StartupReplay
{
public:
  void readParameters(ros::NodeHandle &nh, bool imu_en);
  void initializeLog(const std::string &path);
  bool enabled() const { return enabled_; }
  bool pending() const { return enabled_ && !replay_started_; }
  void cacheImu(const sensor_msgs::Imu::ConstPtr &imu);
  void cacheMeasurement(const LidarMeasureGroup &group);
  void replay(LIVMapper &mapper);
  void finish();

private:
  bool enabled_ = false;
  // Set before executing cached updates: they must publish and must not recache.
  bool replay_started_ = false;
  size_t bytes_ = 0, limit_ = 0;
  std::deque<LidarMeasureGroup> measurements_;
  std::deque<sensor_msgs::Imu::ConstPtr> imu_samples_;
  std::ofstream log_;
};
#endif
