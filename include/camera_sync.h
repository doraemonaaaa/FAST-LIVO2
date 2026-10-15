#ifndef FAST_LIVO_CAMERA_SYNC_H
#define FAST_LIVO_CAMERA_SYNC_H
#include <algorithm>
#include <cmath>
#include <vector>

// Missing views close on sensor-time progress or explicit EOF, never on a
// wall-clock timeout. A stream must be monotonic and within configured lateness.
inline bool cameraBundleReady(double timestamp, const std::vector<double> &first,
    const std::vector<double> &last, double tolerance, double lateness, bool eof)
{
  double watermark = timestamp;
  for (double time : last) watermark = std::max(watermark, time);
  for (size_t k = 0; k < first.size(); ++k)
    if (std::abs(first[k]-timestamp)>tolerance && last[k]<=timestamp+tolerance &&
        watermark<=timestamp+lateness && !eof) return false;
  return true;
}
#endif
