#pragma once
#include "common_lib.h"
#include <memory>

namespace initialization {
// Offline calibration input, independent of the online ESIKF and dataset format.
struct ImuSample {
  using ConstPtr=std::shared_ptr<const ImuSample>;
  double timestamp=0;
  V3D angular_velocity=V3D::Zero(),acceleration=V3D::Zero();
  Eigen::Quaterniond orientation=Eigen::Quaterniond::Identity();
  bool orientation_valid=false;
};
struct TimedPoint {PointType point;double time;};
struct LivoScan {double begin=0,end=0;std::vector<TimedPoint> points;};
inline LivoScan CalibrationScan(const PointCloudXYZI& points,double begin) {
  LivoScan scan;scan.begin=begin;scan.end=begin;scan.points.reserve(points.size());
  for(const auto& point:points) {
    const double time=begin+point.curvature*.001;
    if(!std::isfinite(time) || !point.getVector3fMap().allFinite())
      throw std::runtime_error("Non-finite calibration point");
    scan.end=std::max(scan.end,time);scan.points.push_back({point,time});
  }
  return scan;
}
inline ImuSample::ConstPtr CalibrationImu(const sensor_msgs::Imu& sample,double scale=1.) {
  auto out=std::make_shared<ImuSample>();out->timestamp=sample.header.stamp.toSec();
  out->angular_velocity=V3D(sample.angular_velocity.x,sample.angular_velocity.y,sample.angular_velocity.z);
  out->acceleration=scale*V3D(sample.linear_acceleration.x,sample.linear_acceleration.y,sample.linear_acceleration.z);
  out->orientation=Eigen::Quaterniond(sample.orientation.w,sample.orientation.x,sample.orientation.y,sample.orientation.z);
  out->orientation_valid=sample.orientation_covariance[0]>=0 && out->orientation.coeffs().allFinite() && std::abs(out->orientation.norm()-1)<.01;
  return out;
}
} // namespace initialization
