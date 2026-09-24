#pragma once
#include <cmath>

namespace flivo {
// Equal independent raw-image noise gives residual variance proportional to
// current_gain^2 + reference_gain^2. Divide by two to preserve the existing
// img_point_cov/outlier_threshold calibration when both gains are one.
inline double photometricSqrtInformation(double current_gain, double reference_gain)
{
  if (!std::isfinite(current_gain) || !std::isfinite(reference_gain) ||
      current_gain <= 0.0 || reference_gain <= 0.0) return 0.0;
  const double scale = std::hypot(current_gain, reference_gain);
  const double weight = std::sqrt(2.0) / scale;
  return std::isfinite(weight) ? weight : 0.0;
}
} // namespace flivo
