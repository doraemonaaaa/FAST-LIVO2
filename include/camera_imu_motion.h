#pragma once

#include "common_lib.h"

namespace flivo
{
using CameraPoseJacobian = Eigen::Matrix<double, 6, DIM_STATE>;

class CameraImuMotion
{
public:
  CameraImuMotion() = default;
  CameraImuMotion(const deque<sensor_msgs::Imu::ConstPtr> &samples,
                  double reference_time, double exposure_time, double acceleration_scale);
  bool active() const { return !steps_.empty(); }
  StatesGroup evaluate(const StatesGroup &reference, CameraPoseJacobian *jacobian = nullptr) const;

private:
  struct Step
  {
    double dt;
    V3D angular_velocity;
    V3D acceleration;
  };
  vector<Step, Eigen::aligned_allocator<Step>> steps_;
};
}
