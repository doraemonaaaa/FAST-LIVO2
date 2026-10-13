/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include "IMU_Processing.h"
#include <stdexcept>

void ImuProcess::rewindForStartupReplay(StatesGroup &state, double end, double begin,
    const deque<sensor_msgs::Imu::ConstPtr> &samples)
{
  if (imu_need_init || samples.empty() || !std::isfinite(begin) || end < begin)
    throw std::runtime_error("Startup replay requires initialized state and ordered IMU coverage");
  auto integration = samples;
  if (integration.front()->header.stamp.toSec() > begin)
  {
    sensor_msgs::Imu::Ptr held(new sensor_msgs::Imu(*integration.front()));
    held->header.stamp.fromSec(begin);
    integration.push_front(held);
  }
  // Reverse midpoint/last-sample-hold inertial dynamics over cached IMU times.
  // Add process uncertainty while going backwards; do not treat the terminal
  // filtered covariance as an exact initial-state observation.
  double cursor = end;
  for (int i = int(integration.size()) - 1; i >= 0 && cursor > begin; --i)
  {
    const auto &head = integration[i];
    double left = std::max(begin, head->header.stamp.toSec());
    if (i == 0) left = begin; // At most the missing leading sample interval.
    if (left >= cursor) continue;
    const auto &tail = i + 1 < int(integration.size()) ? integration[i + 1] : head;
    double dt = cursor - left;
    if (dt > .1)
      ROS_WARN("[startup replay] IMU interval %.6f s (%.9f..%.9f): midpoint/hold approximation with process uncertainty",
          dt, left, cursor);
    V3D omega = .5 * V3D(head->angular_velocity.x + tail->angular_velocity.x,
        head->angular_velocity.y + tail->angular_velocity.y,
        head->angular_velocity.z + tail->angular_velocity.z) - state.bias_g;
    V3D acc = .5 * V3D(head->linear_acceleration.x + tail->linear_acceleration.x,
        head->linear_acceleration.y + tail->linear_acceleration.y,
        head->linear_acceleration.z + tail->linear_acceleration.z)
        * G_m_s2 / IMU_mean_acc_norm - state.bias_a;
    M3D rotation = state.rot_end * Exp(omega, -dt);
    V3D acceleration = state.rot_end * acc + state.gravity;
    V3D velocity = state.vel_end - acceleration * dt;
    state.pos_end -= velocity * dt + .5 * acceleration * dt * dt;
    state.vel_end = velocity;
    MD(DIM_STATE, DIM_STATE) F = MD(DIM_STATE, DIM_STATE)::Identity();
    MD(DIM_STATE, DIM_STATE) Q = MD(DIM_STATE, DIM_STATE)::Zero();
    F.block<3,3>(0,0) = Exp(omega, -dt);
    if (ba_bg_est_en) F.block<3,3>(0,10) = -Eye3d * dt;
    F.block<3,3>(3,7) = Eye3d * dt;
    M3D acc_skew; acc_skew << SKEW_SYM_MATRX(acc);
    F.block<3,3>(7,0) = -rotation * acc_skew * dt;
    if (ba_bg_est_en) F.block<3,3>(7,13) = -rotation * dt;
    if (gravity_est_en) F.block<3,3>(7,16) = Eye3d * dt;
    Q.block<3,3>(0,0).diagonal() = cov_gyr * dt * dt;
    Q.block<3,3>(7,7) = rotation * cov_acc.asDiagonal() * rotation.transpose() * dt * dt;
    Q.block<3,3>(10,10).diagonal() = cov_bias_gyr * dt * dt;
    Q.block<3,3>(13,13).diagonal() = cov_bias_acc * dt * dt;
    if (exposure_estimate_en) Q(6,6) = cov_inv_expo * dt * dt;
    MD(DIM_STATE, DIM_STATE) inverse = F.inverse();
    state.cov = inverse * (state.cov + Q) * inverse.transpose();
    state.rot_end = rotation;
    cursor = left;
  }
  if (cursor > begin || !state.pos_end.allFinite() || !state.vel_end.allFinite() || !state.cov.allFinite())
    throw std::runtime_error("Startup replay could not rewind finite state to first input");
  // Change only translation gauge: first input starts at the local origin.
  state.pos_end.setZero();
  pcl_wait_proc.clear(); IMUpose.clear(); cur_pcl_un_->clear();
  last_prop_end_time = begin; imu_time_init = false;
  last_imu = samples.front();
  if (last_imu->header.stamp.toSec() > begin)
  {
    sensor_msgs::Imu::Ptr held(new sensor_msgs::Imu(*last_imu));
    held->header.stamp.fromSec(begin);
    last_imu = held;
  }
  angvel_last = V3D(last_imu->angular_velocity.x, last_imu->angular_velocity.y,
      last_imu->angular_velocity.z) - state.bias_g;
  acc_s_last = state.rot_end * (V3D(last_imu->linear_acceleration.x,
      last_imu->linear_acceleration.y,last_imu->linear_acceleration.z)
      * G_m_s2 / IMU_mean_acc_norm - state.bias_a) + state.gravity;
}

