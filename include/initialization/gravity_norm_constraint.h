#pragma once
#include "common_lib.h"
#include <stdexcept>

namespace livo {
// Soft model constraint, applied once before a formal prediction interval.
// It is not an independent sensor measurement. Update all correlated states,
// use Joseph covariance, and reset right-attitude error after injection.
inline void ConstrainGravityNorm(StatesGroup& state,double sigma,double norm=9.81) {
    const double length=state.gravity.norm();
    if(!std::isfinite(length) || length<1e-6 || !std::isfinite(sigma) || sigma<=0 ||
       !state.cov.allFinite()) throw std::runtime_error("invalid gravity norm constraint state");
    using Matrix=Eigen::MatrixXd;
    using Vector=Eigen::VectorXd;
    const int dimension=state.dimension();
    Eigen::MatrixXd H=Eigen::MatrixXd::Zero(1,dimension);
    H.block<1,3>(0,16)=state.gravity.transpose()/length;
    const Matrix P=(.5*(state.cov+state.cov.transpose())).eval();
    const double variance=(H*P*H.transpose())(0,0)+sigma*sigma;
    if(!std::isfinite(variance) || variance<=0) throw std::runtime_error("invalid gravity norm innovation variance");
    const Vector K=P*H.transpose()/variance;
    const Vector correction=K*(norm-length);
    const Matrix A=Matrix::Identity(dimension,dimension)-K*H;
    Matrix posterior=A*P*A.transpose()+sigma*sigma*K*K.transpose();
    state+=correction;
    const V3D angle=correction.head<3>();const double theta=angle.norm();
    M3D skew;
    skew<<0,-angle.z(),angle.y(),angle.z(),0,-angle.x(),-angle.y(),angle.x(),0;
    const M3D Jr=theta<1e-5 ? (M3D::Identity()-.5*skew+skew*skew/6.).eval()
        : (M3D::Identity()-(1.-std::cos(theta))/(theta*theta)*skew+
           (theta-std::sin(theta))/(theta*theta*theta)*skew*skew).eval();
    Matrix reset=Matrix::Identity(dimension,dimension);reset.block<3,3>(0,0)=Jr;
    posterior=reset*posterior*reset.transpose();
    state.cov=(.5*(posterior+posterior.transpose())).eval();
}
} // namespace livo
