#pragma once
#include <sensor_msgs/Imu.h>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <stdexcept>

namespace imu_support {
constexpr double kMaximumGap=.025;
enum class Coverage {Covered, MissingStart, MissingEnd, Gap, Invalid};
inline const char* Reason(Coverage s) {
  switch(s) {
    case Coverage::Covered:return "covered";
    case Coverage::MissingStart:return "missing_imu_start";
    case Coverage::MissingEnd:return "missing_imu_end";
    case Coverage::Gap:return "imu_gap";
    default:return "invalid_imu_interval";
  }
}
template<class Samples,class Time>
Coverage Check(const Samples& samples,double begin,double end,Time time,
               double maximum_gap=kMaximumGap) {
  if(!std::isfinite(begin) || !std::isfinite(end) || end<begin || maximum_gap<=0)
    return Coverage::Invalid;
  if(samples.empty())return Coverage::MissingEnd;
  if(time(samples.front())>begin)return Coverage::MissingStart;
  if(time(samples.back())<end)return Coverage::MissingEnd;
  auto next=std::upper_bound(samples.begin(),samples.end(),begin,
    [&](double t,const auto& sample){return t<time(sample);});
  if(next==samples.begin())return Coverage::MissingStart;
  auto previous=next-1;
  for(;next!=samples.end() && time(*previous)<end;++previous,++next) {
    const double dt=time(*next)-time(*previous);
    if(!std::isfinite(dt) || dt<=0)return Coverage::Invalid;
    if(dt>maximum_gap)return Coverage::Gap;
  }
  return Coverage::Covered;
}
using Samples=std::deque<sensor_msgs::Imu::ConstPtr>;
inline double Time(const sensor_msgs::Imu::ConstPtr& s){return s->header.stamp.toSec();}
inline Coverage Check(const Samples& s,double begin,double end,double gap=kMaximumGap) {
  return Check(s,begin,end,Time,gap);
}
// Only interpolate BETWEEN acquired samples. Never copy an endpoint outside
// its support or fabricate a constant-value continuation.
inline sensor_msgs::Imu::ConstPtr At(const Samples& samples,double time) {
  if(Check(samples,time,time,std::numeric_limits<double>::infinity())!=Coverage::Covered)
    throw std::runtime_error("IMU interpolation outside measured coverage");
  auto next=std::lower_bound(samples.begin(),samples.end(),time,
    [](const auto& s,double t){return Time(s)<t;});
  if(next!=samples.end() && Time(*next)==time)return *next;
  const auto& a=**(next-1);const auto& b=**next;
  const double alpha=(time-a.header.stamp.toSec())/(b.header.stamp.toSec()-a.header.stamp.toSec());
  sensor_msgs::Imu::Ptr out(new sensor_msgs::Imu(a));out->header.stamp.fromSec(time);
  for(int j=0;j<3;++j) {
    double* w[]={&out->angular_velocity.x,&out->angular_velocity.y,&out->angular_velocity.z};
    double* v[]={&out->linear_acceleration.x,&out->linear_acceleration.y,&out->linear_acceleration.z};
    const double wa[]={a.angular_velocity.x,a.angular_velocity.y,a.angular_velocity.z};
    const double wb[]={b.angular_velocity.x,b.angular_velocity.y,b.angular_velocity.z};
    const double va[]={a.linear_acceleration.x,a.linear_acceleration.y,a.linear_acceleration.z};
    const double vb[]={b.linear_acceleration.x,b.linear_acceleration.y,b.linear_acceleration.z};
    *w[j]=(1-alpha)*wa[j]+alpha*wb[j];*v[j]=(1-alpha)*va[j]+alpha*vb[j];
  }
  Eigen::Quaterniond qa(a.orientation.w,a.orientation.x,a.orientation.y,a.orientation.z);
  Eigen::Quaterniond qb(b.orientation.w,b.orientation.x,b.orientation.y,b.orientation.z);
  if(a.orientation_covariance[0]>=0 && b.orientation_covariance[0]>=0 &&
     qa.coeffs().allFinite() && qb.coeffs().allFinite() && qa.norm()>1e-6 && qb.norm()>1e-6) {
    const auto q=qa.normalized().slerp(alpha,qb.normalized());
    out->orientation.w=q.w();out->orientation.x=q.x();out->orientation.y=q.y();out->orientation.z=q.z();
  } else out->orientation_covariance[0]=-1;
  return out;
}
inline Samples Clip(const Samples& samples,double begin,double end,
                    double gap=std::numeric_limits<double>::infinity()) {
  const auto support=Check(samples,begin,end,gap);
  if(support!=Coverage::Covered)throw std::runtime_error(Reason(support));
  Samples out;out.push_back(At(samples,begin));
  for(const auto& s:samples)if(Time(s)>begin && Time(s)<end)out.push_back(s);
  if(end>begin)out.push_back(At(samples,end));
  return out;
}
} // namespace imu_support
