#ifndef FLIVO_PATCH_SCORE_H
#define FLIVO_PATCH_SCORE_H
#include <algorithm>
#include <cmath>
#include <numeric>
namespace flivo {
inline double cachedPatchMean(const float *patch, int size, float &cached) {
  if (!std::isfinite(cached) || std::abs(cached) < 1e-6)
    cached = static_cast<float>(std::accumulate(patch, patch + size, 0.0) / size);
  return cached;
}
// Signed, pairwise NCC. A constant patch cannot establish a correspondence.
inline bool patchNCC(const float *a, const float *b, int size,
                     double mean_a, double mean_b, double &ncc) {
  double numerator = 0., var_a = 0., var_b = 0.;
  ncc = 0.;
  for (int i = 0; i < size; ++i) {
    const double da = a[i] - mean_a, db = b[i] - mean_b;
    numerator += da * db; var_a += da * da; var_b += db * db;
  }
  if (!std::isfinite(numerator) || !std::isfinite(var_a) || !std::isfinite(var_b) ||
      var_a <= 1e-12 || var_b <= 1e-12) return false;
  ncc = std::max(-1., std::min(1., numerator / std::sqrt(var_a * var_b)));
  return std::isfinite(ncc);
}
}
#endif
