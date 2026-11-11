#include "initialization/gravity_norm_constraint.h"
#include <Eigen/Cholesky>
#include <iostream>
#include <stdexcept>
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
int main(){
    StatesGroup state;state.configureCameras(4);state.gravity=V3D(.1,-.2,-9.73);
    state.cov.setIdentity();state.cov*=.01;
    state.cov(15,18)=state.cov(18,15)=.005;
    state.cov(0,18)=state.cov(18,0)=.001;
    state.cov(19,18)=state.cov(18,19)=.001;
    const double exposure_before=state.extra_exposures[0];
    const double old=std::abs(state.gravity.norm()-9.81);
    livo::ConstrainGravityNorm(state,.02);
    require(std::abs(state.gravity.norm()-9.81)<old,"norm error did not decrease");
    require(std::abs(state.extra_exposures[0]-exposure_before)>1e-5,"correlated camera exposure was not updated");
    require(state.bias_a.norm()>1e-5,"correlated bias was not updated");
    require(state.cov.isApprox(state.cov.transpose(),1e-12),"covariance not symmetric");
    require(Eigen::LLT<decltype(state.cov)>(state.cov).info()==Eigen::Success,"covariance not positive definite");
    require((state.rot_end.transpose()*state.rot_end-M3D::Identity()).norm()<1e-12,"attitude reset invalid");
    StatesGroup uncorrelated;uncorrelated.gravity=V3D(0,0,-9.7);uncorrelated.cov.setIdentity();
    livo::ConstrainGravityNorm(uncorrelated,.02);
    require(uncorrelated.bias_a.norm()==0,"uncorrelated bias changed");
    std::cout<<"PASS gravity model correction, correlated states, Joseph covariance and attitude reset\n";
}
