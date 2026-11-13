#pragma once

#include <Eigen/Dense>
#include <opencv2/core.hpp>
#include <cmath>

namespace flivo
{
struct GraySample
{
  double value, dx, dy;
};

// Coordinates and derivatives are in original-image pixels. A coarse patch
// changes sample spacing, never the four neighbours of an individual sample.
// Caller checks that (floor(x)+1, floor(y)+1) is inside the grayscale image.
inline GraySample sampleGray(const cv::Mat &image, double x, double y)
{
  const int ix = static_cast<int>(std::floor(x)), iy = static_cast<int>(std::floor(y));
  const double u = x - ix, v = y - iy;
  const auto *top = image.ptr<unsigned char>(iy) + ix;
  const auto *bottom = image.ptr<unsigned char>(iy + 1) + ix;
  return {(1-u)*(1-v)*top[0] + u*(1-v)*top[1] + (1-u)*v*bottom[0] + u*v*bottom[1],
          (1-v)*(top[1]-top[0]) + v*(bottom[1]-bottom[0]),
          (1-u)*(bottom[0]-top[0]) + u*(bottom[1]-top[1])};
}

// Eliminate one unconstrained additive brightness offset per patch. Project
// BOTH residuals and Jacobians onto the zero-mean subspace: centering only r
// would give a different objective and an incorrect Kalman update.
inline double centerPatchSystem(Eigen::MatrixXd &H, Eigen::VectorXd &r, int first, int count)
{
  auto rows = H.middleRows(first, count);
  const Eigen::RowVectorXd mean_H = rows.colwise().mean();
  rows.rowwise() -= mean_H;
  auto residuals = r.segment(first, count);
  const double mean_r = residuals.mean();
  residuals.array() -= mean_r;
  return residuals.squaredNorm();
}

inline double centeredPatchError(const float *reference, const float *current, int count,
                                 double reference_gain = 1, double current_gain = 1, double weight = 1)
{
  double mean = 0;
  for (int i = 0; i < count; ++i) mean += current_gain*current[i] - reference_gain*reference[i];
  mean /= count;
  double error = 0;
  for (int i = 0; i < count; ++i)
  {
    const double r = weight*(current_gain*current[i] - reference_gain*reference[i] - mean);
    error += r*r;
  }
  return error;
}
} // namespace flivo
