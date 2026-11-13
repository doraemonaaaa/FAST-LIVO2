#pragma once
#include "initialization/imu_time_offset/imu_time_offset.hpp"
#include <ros/node_handle.h>
#include <fstream>
class LIVMapper;
// Stage 2 only. Stage 1 is the independent offline batch executable.
class StagedInitialization {
public:
  void configure(ros::NodeHandle& nh);
  bool advance(LIVMapper& mapper);
  void observe(const LidarMeasureGroup& group);
  void observeScan(const PointCloudXYZI& points,double begin);
private:
  std::ofstream states_;
  std::string phase_,output_;
  bool enabled_=false;
  int workers_=initialization::kDefaultInitializationThreads;
  double last_=-1,ready_=-1;
  std::vector<initialization::LivoScan> scans_;
  std::vector<initialization::ImuSample::ConstPtr> imus_;
  std::vector<initialization::OffsetInitializationState> snapshots_;
  double chooseOffset(LIVMapper& mapper);
};
