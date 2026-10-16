#pragma once
#include "common_lib.h"
#include "flivo_trace.h"
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace initialization_checkpoint
{
template<class Matrix> void writeMatrix(std::ostream &out, const Matrix &m)
{
  for (int r=0; r<m.rows(); ++r)
    for (int c=0; c<m.cols(); ++c) out << std::setprecision(17) << m(r,c) << ' ';
  out << '\n';
}
template<class Matrix> void readMatrix(std::istream &in, Matrix &m)
{
  for (int r=0; r<m.rows(); ++r)
    for (int c=0; c<m.cols(); ++c)
      if (!(in >> m(r,c)) || !std::isfinite(m(r,c)))
        throw std::runtime_error("Invalid initialization checkpoint matrix");
}
inline void writeState(std::ostream &out, const StatesGroup &s)
{
  out << s.dimension() << '\n';
  writeMatrix(out,s.rot_end); writeMatrix(out,s.pos_end); writeMatrix(out,s.vel_end);
  writeMatrix(out,s.bias_g); writeMatrix(out,s.bias_a); writeMatrix(out,s.gravity);
  out << std::setprecision(17) << s.inv_expo_time << '\n';
  writeMatrix(out,s.extra_exposures); writeMatrix(out,s.cov);
}
inline void readState(std::istream &in, StatesGroup &s)
{
  int dimension=0;
  if (!(in >> dimension) || dimension!=s.dimension())
    throw std::runtime_error("Initialization checkpoint camera/state dimension mismatch");
  readMatrix(in,s.rot_end); readMatrix(in,s.pos_end); readMatrix(in,s.vel_end);
  readMatrix(in,s.bias_g); readMatrix(in,s.bias_a); readMatrix(in,s.gravity);
  if (!(in >> s.inv_expo_time) || !std::isfinite(s.inv_expo_time))
    throw std::runtime_error("Invalid initialization checkpoint exposure");
  readMatrix(in,s.extra_exposures); readMatrix(in,s.cov);
}
inline void writeImu(std::ostream &out, const sensor_msgs::Imu &m)
{
  out << m.header.seq << ' ' << m.header.stamp.toNSec() << ' ' << std::quoted(m.header.frame_id) << '\n';
  for (double x : {m.orientation.x,m.orientation.y,m.orientation.z,m.orientation.w,
       m.angular_velocity.x,m.angular_velocity.y,m.angular_velocity.z,
       m.linear_acceleration.x,m.linear_acceleration.y,m.linear_acceleration.z}) out << std::setprecision(17) << x << ' ';
  for (auto *a : {&m.orientation_covariance,&m.angular_velocity_covariance,&m.linear_acceleration_covariance})
    for (double x : *a) out << std::setprecision(17) << x << ' ';
  out << '\n';
}
inline sensor_msgs::Imu::Ptr readImu(std::istream &in)
{
  sensor_msgs::Imu::Ptr m(new sensor_msgs::Imu()); uint64_t stamp;
  if (!(in >> m->header.seq >> stamp >> std::quoted(m->header.frame_id)))
    throw std::runtime_error("Invalid initialization checkpoint IMU header");
  m->header.stamp.fromNSec(stamp);
  for (double *x : {&m->orientation.x,&m->orientation.y,&m->orientation.z,&m->orientation.w,
       &m->angular_velocity.x,&m->angular_velocity.y,&m->angular_velocity.z,
       &m->linear_acceleration.x,&m->linear_acceleration.y,&m->linear_acceleration.z})
    if (!(in >> *x) || !std::isfinite(*x)) throw std::runtime_error("Invalid initialization checkpoint IMU");
  for (auto *a : {&m->orientation_covariance,&m->angular_velocity_covariance,&m->linear_acceleration_covariance})
    for (double &x : *a) if (!(in >> x) || !std::isfinite(x)) throw std::runtime_error("Invalid IMU covariance");
  return m;
}
inline uint64_t hash(const std::string &value) { return flivo_trace::hashBytes(value.data(),value.size()); }
}
