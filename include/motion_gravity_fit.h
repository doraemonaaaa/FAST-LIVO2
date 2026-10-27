#pragma once
#include "common_lib.h"
#include <Eigen/Eigenvalues>

namespace motion_initialization {
using Matrix6 = Eigen::Matrix<double,6,6>;
using Vector6 = Eigen::Matrix<double,6,1>;
using Matrix5 = Eigen::Matrix<double,5,5>;
using Vector5 = Eigen::Matrix<double,5,1>;
struct GravityFit {
  bool valid=false;
  V3D initial=V3D::Zero(), gravity=V3D::Zero();
  // Endpoint velocity and gravity, including their cross covariance.
  Matrix6 covariance=Matrix6::Identity();
  double error=0;
  int inliers=0;
};

template<int N>
bool solve(const Eigen::Matrix<double,N,N> &H, const Eigen::Matrix<double,N,1> &b,
           Eigen::Matrix<double,N,1> &x) {
  if(!H.allFinite() || !b.allFinite() || (H.diagonal().array()<=1e-12).any())return false;
  const Eigen::Matrix<double,N,N> scale=H.diagonal().cwiseSqrt().cwiseInverse().asDiagonal();
  const Eigen::Matrix<double,N,N> normalized=scale*H*scale;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,N,N>> es(normalized);
  if(es.info()!=Eigen::Success || es.eigenvalues().minCoeff()<1e-6*es.eigenvalues().maxCoeff())return false;
  x=scale*normalized.ldlt().solve(scale*b);
  return x.allFinite();
}

inline Eigen::Matrix<double,3,2> tangent(const V3D &gravity) {
  const V3D direction=gravity.normalized();
  const V3D axis=std::abs(direction.z())<.9?V3D::UnitZ():V3D::UnitX();
  Eigen::Matrix<double,3,2> B;
  B.col(0)=direction.cross(axis).normalized();B.col(1)=direction.cross(B.col(0));
  return B;
}
inline double weight(double whitened_norm) {
  return std::min(1.,2.5/std::max(1e-12,whitened_norm));
}
inline double huber(double whitened_norm) {
  return whitened_norm<=2.5?.5*whitened_norm*whitened_norm:2.5*(whitened_norm-1.25);
}

// dp already includes gravity0. Estimate p-dp = t*v0 + .5*t^2*(g-gravity0).
// Rotations remain those determined by the independent LiDAR pose bootstrap.
template<class Observations>
GravityFit fitGravity(const Observations &obs, const V3D &gravity0, double end_time) {
  GravityFit out;
  if(obs.size()<5 || !gravity0.allFinite() || gravity0.norm()<1 || !std::isfinite(end_time))return out;
  double min_time=1e100,max_time=0;
  for(const auto &o:obs) {
    if(!std::isfinite(o.dt) || o.dt<=0 || !o.position.allFinite() || !o.dp.allFinite() || !o.information.allFinite())return out;
    Eigen::SelfAdjointEigenSolver<M3D> es(o.information);
    if(es.info()!=Eigen::Success || es.eigenvalues().minCoeff()<=0)return out;
    min_time=std::min(min_time,o.dt);max_time=std::max(max_time,o.dt);
  }
  if(max_time<.8-1e-9 || max_time-min_time<.5)return out;
  Vector6 linear=Vector6::Zero();
  for(int iteration=0;iteration<10;++iteration) {
    Matrix6 H=Matrix6::Zero();Vector6 b=Vector6::Zero();
    for(const auto &o:obs) {
      Eigen::Matrix<double,3,6> A;
      A.leftCols<3>()=M3D::Identity()*o.dt;
      A.rightCols<3>()=M3D::Identity()*(.5*o.dt*o.dt);
      const V3D y=o.position-o.dp, r=y-A*linear;
      const double w=iteration?weight(std::sqrt(std::max(0.,r.dot(o.information*r)))):1.;
      H+=A.transpose()*(w*o.information)*A;b+=A.transpose()*(w*o.information)*y;
    }
    Vector6 next;if(!solve<6>(H,b,next))return out;
    const double change=(next-linear).norm();linear=next;if(change<1e-9)break;
  }
  out.initial=linear.head<3>();out.gravity=gravity0+linear.tail<3>();
  if(out.gravity.norm()<1 || !out.gravity.allFinite())return out;
  out.gravity*=9.81/out.gravity.norm();
  auto cost=[&](const V3D &v,const V3D &g) {
    double total=0;
    for(const auto &o:obs) {
      const V3D r=o.position-o.dp-v*o.dt-.5*(g-gravity0)*o.dt*o.dt;
      total+=huber(std::sqrt(std::max(0.,r.dot(o.information*r))));
    }
    return total;
  };
  auto equations=[&](Matrix5 &H,Vector5 &b) {
    H.setZero();b.setZero();const auto B=tangent(out.gravity);
    for(const auto &o:obs) {
      Eigen::Matrix<double,3,5> A;
      A.leftCols<3>()=M3D::Identity()*o.dt;A.rightCols<2>()=.5*o.dt*o.dt*B;
      const V3D r=o.position-o.dp-out.initial*o.dt-.5*(out.gravity-gravity0)*o.dt*o.dt;
      const M3D W=weight(std::sqrt(std::max(0.,r.dot(o.information*r))))*o.information;
      H+=A.transpose()*W*A;b+=A.transpose()*W*r;
    }
  };
  for(int iteration=0;iteration<20;++iteration) {
    Matrix5 H;Vector5 b,step;equations(H,b);if(!solve<5>(H,b,step))return out;
    if(step.norm()<1e-9)break;
    const auto B=tangent(out.gravity);const double previous=cost(out.initial,out.gravity);
    bool accepted=false;
    for(double alpha=1;alpha>=1./128;alpha*=.5) {
      const V3D v=out.initial+alpha*step.head<3>();
      V3D g=out.gravity+alpha*B*step.tail<2>();g*=9.81/g.norm();
      if(cost(v,g)<=previous) {out.initial=v;out.gravity=g;accepted=true;break;}
    }
    if(!accepted)break;
  }
  double squared=0,chi_squared=0;
  for(const auto &o:obs) {
    const V3D r=o.position-o.dp-out.initial*o.dt-.5*(out.gravity-gravity0)*o.dt*o.dt;
    const double whitened=std::sqrt(std::max(0.,r.dot(o.information*r)));
    if(whitened<=2.5) {++out.inliers;squared+=r.squaredNorm();chi_squared+=whitened*whitened;}
  }
  if(out.inliers<5 || out.inliers<.7*obs.size())return out;
  out.error=std::sqrt(squared/out.inliers);
  Matrix5 H;Vector5 b,unused;equations(H,b);if(!solve<5>(H,b,unused))return out;
  const Matrix5 covariance=H.ldlt().solve(Matrix5::Identity())*
      std::max(1.,chi_squared/std::max(1,3*out.inliers-5));
  Eigen::Matrix<double,6,5> lift=Eigen::Matrix<double,6,5>::Zero();
  lift.topLeftCorner<3,3>().setIdentity();
  lift.topRightCorner<3,2>()=end_time*tangent(out.gravity);
  lift.bottomRightCorner<3,2>()=tangent(out.gravity);
  out.covariance=lift*covariance*lift.transpose();
  // Conservative velocity floor; a nonzero gravity floor keeps the KF's
  // full Cartesian covariance invertible despite the fixed gravity norm.
  out.covariance.topLeftCorner<3,3>()+=M3D::Identity()*.01;
  out.covariance.bottomRightCorner<3,3>()+=M3D::Identity()*.0001;
  out.valid=out.initial.allFinite() && out.gravity.allFinite() && out.covariance.allFinite() && out.error<.08;
  return out;
}
} // namespace motion_initialization
