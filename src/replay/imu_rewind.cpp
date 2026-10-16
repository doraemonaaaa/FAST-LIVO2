/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include "IMU_Processing.h"
#include "state_numerics.h"
#include "replay/initialization_checkpoint.h"
#include <stdexcept>

void ImuProcess::writeInitialization(std::ostream &out) const
{
  using namespace initialization_checkpoint;
  if (imu_need_init || !last_imu) throw std::runtime_error("Cannot checkpoint uninitialized IMU");
  out << std::setprecision(17) << IMU_mean_acc_norm << ' ' << cov_inv_expo << ' '
      << last_prop_end_time << ' ' << first_lidar_time << ' ' << imu_time_init << '\n';
  for (const auto *v : {&cov_acc,&cov_gyr,&cov_bias_acc,&cov_bias_gyr,
                        &mean_acc,&mean_gyr,&angvel_last,&acc_s_last}) writeMatrix(out,*v);
  writeImu(out,*last_imu);
}

void ImuProcess::readInitialization(std::istream &in)
{
  using namespace initialization_checkpoint;
  if (!(in >> IMU_mean_acc_norm >> cov_inv_expo >> last_prop_end_time >> first_lidar_time >> imu_time_init)
      || !std::isfinite(IMU_mean_acc_norm) || IMU_mean_acc_norm<=0 || !std::isfinite(cov_inv_expo)
      || !std::isfinite(last_prop_end_time) || !std::isfinite(first_lidar_time))
    throw std::runtime_error("Invalid initialization checkpoint IMU parameters");
  for (auto *v : {&cov_acc,&cov_gyr,&cov_bias_acc,&cov_bias_gyr,
                  &mean_acc,&mean_gyr,&angvel_last,&acc_s_last}) readMatrix(in,*v);
  last_imu=readImu(in);
  pcl_wait_proc.clear(); IMUpose.clear(); cur_pcl_un_->clear();
  imu_need_init=false;
}

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
    MatrixXd F = MatrixXd::Identity(state.dimension(), state.dimension());
    MatrixXd Q = MatrixXd::Zero(state.dimension(), state.dimension());
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
    if (exposure_estimate_en)
      for (size_t k = 0; k <= size_t(state.extra_exposures.size()); ++k)
        Q(StatesGroup::exposureIndex(k), StatesGroup::exposureIndex(k)) = cov_inv_expo * dt * dt;
    const MD(DIM_STATE,DIM_STATE) base = F.topLeftCorner<DIM_STATE,DIM_STATE>();
    const MD(DIM_STATE,DIM_STATE) inverse = base.inverse();
    const MD(DIM_STATE,DIM_STATE) covariance = state.cov.topLeftCorner<DIM_STATE,DIM_STATE>();
    const MD(DIM_STATE,DIM_STATE) noise = Q.topLeftCorner<DIM_STATE,DIM_STATE>();
    state.cov.topLeftCorner<DIM_STATE,DIM_STATE>() = inverse*(covariance+noise)*inverse.transpose();
    const int extra = state.dimension()-DIM_STATE;
    if (extra)
    {
      state.cov.topRightCorner(DIM_STATE,extra) = (inverse*state.cov.topRightCorner(DIM_STATE,extra)).eval();
      state.cov.bottomLeftCorner(extra,DIM_STATE) = (state.cov.bottomLeftCorner(extra,DIM_STATE)*inverse.transpose()).eval();
      state.cov.bottomRightCorner(extra,extra) += Q.bottomRightCorner(extra,extra);
    }
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

