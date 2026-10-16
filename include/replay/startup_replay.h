#ifndef FAST_LIVO_STARTUP_REPLAY_H
#define FAST_LIVO_STARTUP_REPLAY_H
#include "common_lib.h"
#include <fstream>
#include <ros/node_handle.h>

class LIVMapper;

// First process retains IMU and indices only. A fresh second process restores
// the initialization checkpoint and rereads the original sensor messages.
class StartupReplay
{
public:
  void readParameters(ros::NodeHandle &nh, bool imu_en);
  void initializeLog(const std::string &path);
  bool enabled() const { return phase_ != Phase::Disabled; }
  bool pending() const { return phase_ == Phase::Preparing; }
  bool replaying() const { return phase_ == Phase::Replaying; }
  bool prepared() const { return phase_ == Phase::Prepared; }
  void retainImuForInitialization(const sensor_msgs::Imu::ConstPtr &imu);
  void observeMeasurement(const LidarMeasureGroup &group);
  void prepareNextInput(LidarMeasureGroup &group) const;
  void afterUpdate(LIVMapper &mapper);
  void restore(LIVMapper &mapper);
  void finish();
private:
  struct Index { double begin,end; uint64_t fingerprint; size_t points,images; int flag; };
  Index index(const LidarMeasureGroup &group) const;
  enum class Phase { Disabled, Preparing, Prepared, Replaying, Complete };
  Phase phase_=Phase::Disabled;
  std::string checkpoint_, identity_;
  size_t bytes_=0, limit_=0, verified_=0;
  double begin_=-1, end_=0;
  std::vector<Index> indices_;
  std::deque<sensor_msgs::Imu::ConstPtr> imu_samples_;
  std::ofstream log_;
  std::ofstream input_log_;
};
#endif
