#pragma once

#include <algorithm>
#include <cmath>
#include <Eigen/Core>
#include <vikit/abstract_camera.h>

namespace flivo {

// This frontend uses forward-facing bearings. In particular, the vikit
// equidistant implementation divides by Z and folds negative-Z rays across
// the image. Do not sample that branch or differentiate across Z=0.
inline bool validProjectionRay(const Eigen::Vector3d& p)
{
  const double range = p.norm();
  return p.allFinite() && std::isfinite(range) && range > 1e-12 &&
         p.z() > 1e-8 * range;
}

// Differentiate the very same virtual projection used for pixel sampling.
// This includes config-selected distortion, focal lengths and image scaling;
// no second camera-model switch or duplicate calibration is maintained here.
// Normalize the ray to make the step independent of scene units, then apply
// the 1/range chain factor. Six projections are needed per point, not per pixel.
inline bool projectionJacobian(const vk::AbstractCamera& camera,
                               const Eigen::Vector3d& p,
                               Eigen::Matrix<double, 2, 3>& jacobian)
{
  jacobian.setZero();
  if (!validProjectionRay(p)) return false;
  const double range = p.norm();
  const Eigen::Vector3d ray = p / range;
  for (int axis = 0; axis < 3; ++axis)
  {
    const double step = axis == 2 ? std::min(1e-6, 0.25 * ray.z()) : 1e-6;
    Eigen::Vector3d plus = ray, minus = ray;
    plus[axis] += step;
    minus[axis] -= step;
    const Eigen::Vector2d upper = camera.world2cam(plus);
    const Eigen::Vector2d lower = camera.world2cam(minus);
    if (!upper.allFinite() || !lower.allFinite())
    {
      jacobian.setZero();
      return false;
    }
    jacobian.col(axis) = (upper - lower) / (2.0 * step * range);
  }
  if (!jacobian.allFinite())
  {
    jacobian.setZero();
    return false;
  }
  return true;
}

} // namespace flivo
