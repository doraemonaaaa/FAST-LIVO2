#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

// Corrected acquisition order; tie breaking follows configured camera order.
template<class Images>
std::vector<size_t> cameraEventOrder(const Images &images) {
  std::vector<size_t> order(images.size());std::iota(order.begin(),order.end(),0);
  for(const auto &image:images)
    if(!std::isfinite(image.timestamp))throw std::runtime_error("Non-finite camera event time");
  std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b) {
    if(images[a].timestamp!=images[b].timestamp)return images[a].timestamp<images[b].timestamp;
    return images[a].camera_index<images[b].camera_index;
  });
  return order;
}

// An empty stream can still deliver an earlier image. Wait for its progress,
// or close that wait using bounded sensor-time lateness / explicit EOF.
inline bool cameraEventReady(double timestamp,const std::vector<double> &first,
    const std::vector<double> &last,double sensor_time,double lateness,bool eof) {
  for(size_t k=0;k<first.size();++k)
    if(!std::isfinite(first[k]) && last[k]<timestamp && !eof && sensor_time<timestamp+lateness)
      return false;
  return true;
}
