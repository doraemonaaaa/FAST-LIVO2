#pragma once

// Pending point times are relative to the previous image boundary. A LiDAR
// scan can span several images, so consuming its entire remainder is unsafe.
template <typename PointRange>
void splitPendingPointsAtImageTime(PointRange &pending, PointRange &current,
                                  double interval_ms)
{
  PointRange carry;
  carry.swap(pending);
  current.clear();
  current.reserve(carry.size());
  pending.reserve(carry.size());
  for (auto point : carry)
  {
    if (point.curvature < interval_ms)
      current.push_back(point);
    else
    {
      point.curvature -= interval_ms;
      pending.push_back(point);
    }
  }
}
