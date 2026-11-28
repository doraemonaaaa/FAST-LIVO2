#include "camera_imu_motion.h"
#include "imu_support.h"

namespace flivo
{
namespace
{
M3D skew(const V3D &v)
{
  M3D result;
  result << SKEW_SYM_MATRX(v);
  return result;
}

M3D rightJacobian(const V3D &rotation)
{
  const double angle2 = rotation.squaredNorm();
  const M3D cross = skew(rotation);
  if (angle2 < 1e-8)
    return M3D::Identity() - .5 * cross + (1. / 6.) * cross * cross;
  const double angle = std::sqrt(angle2);
  return M3D::Identity() - ((1. - std::cos(angle)) / angle2) * cross +
      ((angle - std::sin(angle)) / (angle2 * angle)) * cross * cross;
}
}

CameraImuMotion::CameraImuMotion(const imu_support::Samples &samples,
    double reference_time, double exposure_time, double acceleration_scale)
{
  if (reference_time == exposure_time) return;
  const auto interval = imu_support::Clip(samples, std::min(reference_time, exposure_time),
      std::max(reference_time, exposure_time), imu_support::kMaximumGap);
  const bool forward = exposure_time > reference_time;
  for (size_t i = 0; i + 1 < interval.size(); ++i)
  {
    const size_t index = forward ? i : interval.size() - 1 - i;
    const auto &a = *interval[index];
    const auto &b = *interval[forward ? index + 1 : index - 1];
    const V3D w(a.angular_velocity.x + b.angular_velocity.x,
                a.angular_velocity.y + b.angular_velocity.y,
                a.angular_velocity.z + b.angular_velocity.z);
    const V3D acc(a.linear_acceleration.x + b.linear_acceleration.x,
                  a.linear_acceleration.y + b.linear_acceleration.y,
                  a.linear_acceleration.z + b.linear_acceleration.z);
    steps_.push_back({b.header.stamp.toSec() - a.header.stamp.toSec(),
                      .5 * w, (.5 * acceleration_scale) * acc});
  }
}

StatesGroup CameraImuMotion::evaluate(const StatesGroup &reference, CameraPoseJacobian *jacobian) const
{
  StatesGroup pose = reference;
  using Derivative = Eigen::Matrix<double, 3, DIM_STATE>;
  Derivative dr, dp, dv;
  if (jacobian)
  {
    dr.setZero(); dp.setZero(); dv.setZero();
    dr.block<3, 3>(0, 0).setIdentity();
    dp.block<3, 3>(0, 3).setIdentity();
    dv.block<3, 3>(0, 7).setIdentity();
  }
  for (const auto &step : steps_)
  {
    const double dt = step.dt;
    const V3D rate = step.angular_velocity - reference.bias_g;
    const V3D specific_force = step.acceleration - reference.bias_a;
    const V3D rotation = rate * dt;
    const double angle = rotation.norm();
    M3D increment = M3D::Identity();
    if (angle > 0) increment = Eigen::AngleAxisd(angle, rotation / angle).toRotationMatrix();
    const M3D next_rotation = pose.rot_end * increment;
    // The deskew integrator evaluates acceleration at the later endpoint.
    // Reversing the step uses the same endpoint, giving the inverse motion.
    const M3D &acc_rotation = dt > 0 ? next_rotation : pose.rot_end;
    const V3D acceleration = acc_rotation * specific_force + reference.gravity;
    if (jacobian)
    {
      Derivative next_dr = increment.transpose() * dr;
      next_dr.block<3, 3>(0, 10) -= dt * rightJacobian(rotation);
      Derivative da = -acc_rotation * skew(specific_force) * (dt > 0 ? next_dr : dr);
      da.block<3, 3>(0, 13) -= acc_rotation;
      da.block<3, 3>(0, 16).diagonal().array() += 1.;
      dp += dt * dv + (.5 * dt * dt) * da;
      dv += dt * da;
      dr = next_dr;
    }
    pose.pos_end += dt * pose.vel_end + (.5 * dt * dt) * acceleration;
    pose.vel_end += dt * acceleration;
    pose.rot_end = next_rotation;
  }
  if (jacobian)
  {
    jacobian->topRows<3>() = dr;
    jacobian->bottomRows<3>() = dp;
  }
  return pose;
}
}
