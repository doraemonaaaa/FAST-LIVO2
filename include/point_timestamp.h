#pragma once
#include <sensor_msgs/PointCloud2.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace point_time {
struct Bounds { double first, last; };
// Read the same first-point origin used by preprocessing, before filtering.
inline Bounds bounds(const sensor_msgs::PointCloud2& cloud, double max_duration) {
  if (!std::isfinite(max_duration) || max_duration <= 0)
    throw std::runtime_error("invalid maximum scan duration");
  const sensor_msgs::PointField* field = nullptr;
  for (const auto& f : cloud.fields) if (f.name == "timestamp") field = &f;
  if (!field || field->datatype != sensor_msgs::PointField::FLOAT64 || field->count != 1)
    throw std::runtime_error("point_timestamp requires scalar FLOAT64 timestamp in seconds");
  if (!cloud.width || !cloud.height || cloud.point_step < 8 || field->offset > cloud.point_step - 8 ||
      uint64_t(cloud.row_step) < uint64_t(cloud.width) * cloud.point_step ||
      uint64_t(cloud.data.size()) < uint64_t(cloud.row_step) * cloud.height)
    throw std::runtime_error("invalid point cloud timestamp layout");
  const uint16_t endian = 1;
  const bool host_big = *reinterpret_cast<const unsigned char*>(&endian) == 0;
  Bounds result{0,0};
  for (uint32_t row=0; row<cloud.height; ++row)
    for (uint32_t col=0; col<cloud.width; ++col) {
      unsigned char bytes[8];
      std::memcpy(bytes,cloud.data.data()+size_t(row)*cloud.row_step+size_t(col)*cloud.point_step+field->offset,8);
      if (host_big != cloud.is_bigendian) std::reverse(bytes,bytes+8);
      double stamp; std::memcpy(&stamp,bytes,8);
      if (!std::isfinite(stamp) || stamp<=0) throw std::runtime_error("invalid point timestamp");
      if (row==0 && col==0) result={stamp,stamp};
      if (stamp<result.first || stamp-result.first>max_duration)
        throw std::runtime_error("point timestamps precede first point or exceed scan duration");
      result.last=std::max(result.last,stamp);
    }
  return result;
}
inline double anchor(const Bounds& times,double clock_offset) {
  if (!std::isfinite(clock_offset)) throw std::runtime_error("invalid LiDAR clock offset");
  return times.first+clock_offset;
}
inline double rebased_ms(double offset_ms,double scan_start,double propagation_start) {
  return offset_ms+(scan_start-propagation_start)*1000.0;
}
}
