#pragma once

#include "common_lib.h"
#include <ros/node_handle.h>
#include <memory>

class ImuProcess;

// Export additional scans in the IMU frame of a published main LiDAR pose.
class AuxiliaryLidarExport
{
public:
  AuxiliaryLidarExport(ros::NodeHandle &nh, const std::string &directory,
                       double blind_range, bool save_scans);
  ~AuxiliaryLidarExport();
  void observeImu(const sensor_msgs::Imu::ConstPtr &sample);
  void observePose(double time, const StatesGroup &state);
  void drain(const ImuProcess &imu, double estimator_time, bool finished);

private:
  struct Reference
  {
    double time;
    StatesGroup state;
  };
  struct Stream;

  std::deque<sensor_msgs::Imu::ConstPtr> samples_;
  std::deque<Reference> references_;
  std::vector<std::unique_ptr<Stream>> streams_;
};
