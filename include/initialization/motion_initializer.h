#include "imu_support.h"
#pragma once
#include "common_lib.h"
#include "motion_gravity_fit.h"
#include <map>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <Eigen/Eigenvalues>
#include <deque>

// Short-window, attitude-aided LiDAR pose bootstrap. No map/EKF feedback.
// Register endpoint clouds to consecutive scans; fit motion in the window anchor frame.
class MotionInitializer {
public:
  void setBias(const V3D &bg,const V3D &ba,double scale) {bg_=bg;ba_=ba;scale_=scale;}
  using PoseMatrix = Eigen::Matrix<double,6,6>;
  using PoseVector = Eigen::Matrix<double,6,1>;
  explicit MotionInitializer(bool timeout_fallback=false, double velocity_sigma=1.0, bool estimate_gravity=true)
      : timeout_fallback_(timeout_fallback), fallback_velocity_sigma_(std::max(1.0, velocity_sigma)),
        estimate_gravity_(estimate_gravity) {}
  struct Result {
    bool ready=false, timeout_fallback=false;
    V3D velocity=V3D::Zero(), position=V3D::Zero(), gravity=V3D::Zero();
    M3D rotation=M3D::Identity(), velocity_cov=M3D::Identity();
    PoseMatrix pose_covariance=PoseMatrix::Identity();
    PoseMatrix velocity_gravity_covariance=PoseMatrix::Identity();
    bool gravity_estimated=false;
    double gravity_correction=0;
    int gravity_stable_windows=0;
    double elapsed=0, match_ratio=0, residual=0, eigen_ratio=0, fit_error=0, rotation_correction=0;
  };
  struct Observation { double dt; V3D position, dp, dv; M3D information; };
  struct VelocityFit { bool valid=false; V3D initial=V3D::Zero(); M3D covariance=M3D::Identity(); double error=0; };
  static VelocityFit fit(const std::vector<Observation>& obs) {
    VelocityFit out; M3D H=M3D::Zero(); V3D b=V3D::Zero();
    for(const auto &o:obs) { H+=o.dt*o.dt*o.information; b+=o.dt*o.information*(o.position-o.dp); }
    Eigen::SelfAdjointEigenSolver<M3D> es(H);
    if(obs.size()<5 || es.eigenvalues().minCoeff()<=1e-6) return out;
    out.initial=H.ldlt().solve(b); double e=0;
    for(const auto &o:obs) e+=(o.position-o.dp-out.initial*o.dt).squaredNorm();
    out.error=std::sqrt(e/obs.size());
    // Conservatively retain at least 0.1 m/s uncertainty per axis.
    out.covariance=H.inverse()+M3D::Identity()*.01;
    out.valid=out.initial.allFinite() && out.covariance.allFinite() && out.error<.08;
    return out;
  }
  // Both the information and the rotation difference use anchor-frame axes.
  static constexpr double integration_gyro_sigma=0.003490658503988659; // radians, 0.2 degrees
  static M3D integrationRotationGain(const M3D &lidar_information) {
    return (lidar_information+M3D::Identity()/(integration_gyro_sigma*integration_gyro_sigma))
              .ldlt().solve(lidar_information);
  }
  Result update(const PointCloudXYZI &scan, double begin, double end,
                const std::deque<sensor_msgs::Imu::ConstPtr> &imus,
                const M3D &R_LI, const V3D &t_LI) {
    Result out;
    if(imus.empty() || scan.empty()) { resetWindow(); return out; }
    if(imus.back()->header.stamp.toNSec()<ros::Time().fromSec(end).toNSec())
      throw std::runtime_error("motion initialization would extrapolate beyond acquired IMU");
    const auto &m=imus.back(); Eigen::Quaterniond q(m->orientation.w,m->orientation.x,m->orientation.y,m->orientation.z);
    if(m->orientation_covariance[0]<0 || !q.coeffs().allFinite() || std::abs(q.norm()-1)>.01)
      throw std::runtime_error("motion initialization: invalid IMU attitude");
    q.normalize();
    if(first_attempt<0)first_attempt=end;
    if(start<0) {
      start=end; anchor=q; previous_time=end; last_scan_time=end;
      gyro_rotation_.setIdentity();previous_omega_=V3D(m->angular_velocity.x,m->angular_velocity.y,m->angular_velocity.z)-bg_;
      last_sample={end,M3D::Identity(),scale_*V3D(m->linear_acceleration.x,m->linear_acceleration.y,m->linear_acceleration.z)-ba_};
      if(!last_sample.acceleration.allFinite()) throw std::runtime_error("motion initialization: invalid acceleration");
      previous_acc=last_sample.acceleration+anchor.conjugate()*V3D(0,0,-9.81);
    }
    out.elapsed=end-start;
    if(end-first_attempt>10) {
      if(!timeout_fallback_)
        throw std::runtime_error("motion initialization: no observable consistent velocity within 10 seconds");
      // Restart the local frame at the current scan. Do not use unaccepted motion fits.
      out.ready=true; out.timeout_fallback=true;
      out.gravity=q.conjugate()*V3D(0,0,-9.81);
      out.velocity_cov=M3D::Identity()*fallback_velocity_sigma_*fallback_velocity_sigma_;
      return out;
    }
    V3D gravity=anchor.conjugate()*V3D(0,0,-9.81);
    // Retain only IMU samples since the last accepted pose. Their acceleration
    // is integrated once the LiDAR endpoint rotation for this interval is known.
    for(const auto &im:imus) {
      double time=im->header.stamp.toSec(); if(time<=previous_time) continue;
      double dt=time-previous_time;
      if(dt>imu_support::kMaximumGap) {
        ROS_WARN("[motion init] IMU gap %.6f s: discard startup window and restart",dt);
        resetWindow();
        return update(scan,begin,end,imus,R_LI,t_LI);
      }
      const V3D acceleration=scale_*V3D(im->linear_acceleration.x,im->linear_acceleration.y,im->linear_acceleration.z)-ba_;
      const V3D omega=V3D(im->angular_velocity.x,im->angular_velocity.y,im->angular_velocity.z)-bg_;
      if(!omega.allFinite()) throw std::runtime_error("motion initialization: invalid angular velocity");
      const V3D midpoint=.5*(previous_omega_+omega);
      gyro_rotation_=gyro_rotation_*Exp(midpoint,dt);previous_omega_=omega;
      if(!acceleration.allFinite()) throw std::runtime_error("motion initialization: invalid acceleration");
      if(!omega.allFinite()) throw std::runtime_error("motion initialization: invalid angular velocity");
      pending_imus.push_back({time,gyro_rotation_,acceleration});
      previous_time=time;
    }
    // End-of-scan attitude, rotational deskew, and the latest independent velocity estimate.
    const M3D imu_rotation=gyro_rotation_;
    M3D relative_rotation=last_imu_rotation.transpose()*imu_rotation;
    const M3D endpoint_prior=last_rotation*relative_rotation;
    const V3D endpoint_velocity=endpoint_prior.transpose()*velocity_guess;
    V3D omega=V3D(m->angular_velocity.x,m->angular_velocity.y,m->angular_velocity.z)-bg_;
    PointCloudXYZI::Ptr cloud(new PointCloudXYZI());
    for(const auto &pt:scan) {
      V3D p(pt.x,pt.y,pt.z); if(!p.allFinite() || p.norm()<1 || p.norm()>45) continue;
      double dt=end-(begin+double(pt.curvature)*.001);
      if(dt<-.01 || dt>.3) continue;
      V3D w=-omega*dt; double angle=w.norm();M3D dR=M3D::Identity();if(angle>1e-10)dR=Eigen::AngleAxisd(angle,w/angle).toRotationMatrix();
      p=dR*(R_LI*p+t_LI)-endpoint_velocity*dt;
      PointType v=pt;v.x=p.x();v.y=p.y();v.z=p.z();v.curvature=dt*1000.;cloud->push_back(v);
    }
    // Deterministic voxel centroids; avoid binary PCL/Eigen alignment dependencies.
    struct Cell { V3D sum=V3D::Zero(); int count=0; };
    std::map<std::array<int,3>,Cell> cells;
    for(const auto &p:*cloud) {
      const V3D body(p.x,p.y,p.z);
      std::array<int,3> key{{int(std::floor(body.x()/.25)),int(std::floor(body.y()/.25)),int(std::floor(body.z()/.25))}};
      auto &cell=cells[key];cell.sum+=body;++cell.count;
    }
    PointCloudXYZI::Ptr reduced(new PointCloudXYZI());
    for(const auto &entry:cells) {
      const auto &cell=entry.second;V3D p=cell.sum/cell.count;
      PointType v{};v.x=p.x();v.y=p.y();v.z=p.z();reduced->push_back(v);
    }
    // A rejected pair breaks the measured displacement chain. Keep first_attempt
    // so restarting the local window never extends the overall timeout.
    auto reject_pair = [&]() { resetWindow(); return out; };
    if(reduced->size()<100) return reject_pair();
    if(!target) {
      if(!setTarget(reduced)) return reject_pair();
      return out;
    }
    V3D relative_translation=last_rotation.transpose()*velocity_guess*(end-last_scan_time);
    PoseMatrix H=PoseMatrix::Zero();double rms=0;int count=0;
    size_t plane_candidates=0;
    const M3D gyro_target=last_rotation.transpose()*imu_rotation;
    constexpr double lidar_variance=.03*.03;
    constexpr double registration_sigma=0.05235987755982989; // radians, 3 degrees
    const double prior_weight=lidar_variance/(registration_sigma*registration_sigma);
    PoseMatrix registration_H;
    bool converged=false;
    for(int iteration=0;iteration<=20;++iteration) {
      H.setZero();PoseVector gradient=PoseVector::Zero();double squared=0;count=0;
      plane_candidates=0;
      for(const auto &pt:*reduced) {
        const V3D body(pt.x,pt.y,pt.z);
        const V3D p=relative_rotation*body+relative_translation;
        PointType search{};search.x=p.x();search.y=p.y();search.z=p.z();
        std::vector<int> ids(1);std::vector<float> distances(1);
        if(nearest(search,1,ids,distances)!=1 || distances[0]>.64 || normals[ids[0]].norm()<.5)continue;
        ++plane_candidates;
        const V3D &n=normals[ids[0]];double r=n.dot(p-centers[ids[0]]);if(std::abs(r)>.3)continue;
        double weight=std::min(1.,.05/std::max(1e-9,std::abs(r)));
        // Right rotation increment; translation is in the previous scan frame.
        PoseVector J;J.head<3>()=body.cross(relative_rotation.transpose()*n);J.tail<3>()=n;
        H+=weight*J*J.transpose();gradient+=weight*J*r;squared+=r*r;++count;
      }
      if(count<100)return reject_pair();
      registration_H=H;
      registration_H.topLeftCorner<3,3>().diagonal().array()+=prior_weight;
      rms=std::sqrt(squared/count);
      if(converged || iteration==20) break; // H now describes the returned pose.
      gradient.head<3>()+=prior_weight*Log(M3D(gyro_target.transpose()*relative_rotation));
      PoseVector step=-registration_H.ldlt().solve(gradient);
      if(!step.allFinite() || step.head<3>().norm()>.2 || step.tail<3>().norm()>.75)return reject_pair();
      relative_rotation=relative_rotation*Exp(V3D(step.head<3>()));
      relative_translation+=step.tail<3>();
      converged=step.head<3>().norm()<1e-5 && step.tail<3>().norm()<1e-4;
    }
    out.match_ratio=plane_candidates ? double(count)/plane_candidates : 0.;
    out.residual=rms;
    if(out.match_ratio<.20 || rms>.10)return reject_pair();
    const M3D Hrt=H.topRightCorner<3,3>();
    const M3D Htt=H.bottomRightCorner<3,3>();
    const M3D trans_info=Htt-Hrt.transpose()*registration_H.topLeftCorner<3,3>().ldlt().solve(Hrt);
    Eigen::SelfAdjointEigenSolver<M3D> es(trans_info);
    if(es.info()!=Eigen::Success || es.eigenvalues().maxCoeff()<=0) return reject_pair();
    out.eigen_ratio=es.eigenvalues().minCoeff()/es.eigenvalues().maxCoeff();
    if(out.eigen_ratio<.005 || !setTarget(reduced)) return reject_pair();

    const V3D translation=last_translation+last_rotation*relative_translation;
    const M3D R=(Eigen::Quaterniond(last_rotation)*Eigen::Quaterniond(relative_rotation))
                    .normalized().toRotationMatrix();
    const PoseMatrix pair_covariance=registration_H.ldlt().solve(PoseMatrix::Identity())*(double(count)*.03*.03);
    // Weight each observation by its pairwise registration quality in the anchor axes.
    const M3D info=last_rotation*pair_covariance.bottomRightCorner<3,3>().inverse()*last_rotation.transpose();
    // ESIKF consumes endpoint pose uncertainty, including preceding increments.
    // Compose right attitude errors and additive anchor-frame position errors.
    M3D skew; skew << SKEW_SYM_MATRX(relative_translation);
    PoseMatrix previous_jacobian=PoseMatrix::Identity(), pair_jacobian=PoseMatrix::Identity();
    previous_jacobian.topLeftCorner<3,3>()=relative_rotation.transpose();
    previous_jacobian.bottomLeftCorner<3,3>()=-last_rotation*skew;
    pair_jacobian.bottomRightCorner<3,3>()=last_rotation;
    out.pose_covariance=previous_jacobian*last_pose_covariance*previous_jacobian.transpose()
                        +pair_jacobian*pair_covariance*pair_jacobian.transpose();
    last_pose_covariance=out.pose_covariance;
    // H uses right increments in scan coordinates. Rotate its marginal
    // information into the anchor frame of the left correction R * gyro^-1.
    const M3D marginal=H.topLeftCorner<3,3>()-Hrt*Htt.ldlt().solve(Hrt.transpose());
    const M3D rotation_info=R*(.5*(marginal+marginal.transpose()))*R.transpose()/(count*lidar_variance);
    const V3D difference=Log(M3D(R*imu_rotation.transpose()));
    out.rotation_correction=difference.norm();
    if(out.rotation_correction>0.17453292519943295) return reject_pair();
    const M3D gain=integrationRotationGain(rotation_info);
    const V3D correction=gain*difference;
    if(!correction.allFinite()) return reject_pair();
    const M3D integration_correction=Exp(correction,1.0);
    const M3D fused_rotation=integration_correction*imu_rotation;
    // The ESIKF uses a right rotation error and additive anchor-frame position.
    // Apply the same directional gain to attitude and its position cross-covariance.
    // Keep last_pose_covariance in the raw scan chain for the next registration.
    PoseMatrix to_fused=PoseMatrix::Identity();
    to_fused.topLeftCorner<3,3>()=fused_rotation.transpose()*gain*R;
    const M3D gyro_noise=integration_gyro_sigma*fused_rotation.transpose()*(M3D::Identity()-gain);
    out.pose_covariance=to_fused*last_pose_covariance*to_fused.transpose();
    out.pose_covariance.topLeftCorner<3,3>()+=gyro_noise*gyro_noise.transpose();
    integrateAcceptedMotion(integration_correction,gravity,end);
    last_integration_correction_=integration_correction;
    observations.push_back({out.elapsed,translation,dp,dv,info});
    last_translation=translation;last_scan_time=end;last_rotation=R;last_imu_rotation=imu_rotation;
    out.position=translation;out.rotation=fused_rotation;out.gravity=gravity;
    if(out.elapsed>.15)velocity_guess=(translation-dp)/out.elapsed+dv;
    if(out.elapsed<.8)return out;
    if(estimate_gravity_) {
      const auto fitted=motion_initialization::fitGravity(observations,gravity,out.elapsed);
      out.fit_error=fitted.error;
      if(!fitted.valid) {gravity_stable_windows_=0;return out;}
      const double change=last_fit_gravity_.norm()>0?
          std::acos(std::max(-1.,std::min(1.,last_fit_gravity_.normalized().dot(fitted.gravity.normalized())))):1e100;
      gravity_stable_windows_=change<.017453292519943295?gravity_stable_windows_+1:1;
      last_fit_gravity_=fitted.gravity;
      out.gravity_stable_windows=gravity_stable_windows_;
      out.gravity_correction=std::acos(std::max(-1.,std::min(1.,gravity.normalized().dot(fitted.gravity.normalized()))));
      out.gravity=fitted.gravity;
      out.velocity=fitted.initial+dv+(fitted.gravity-gravity)*out.elapsed;
      out.velocity_cov=fitted.covariance.topLeftCorner<3,3>();
      out.velocity_gravity_covariance=fitted.covariance;
      velocity_guess=out.velocity;
      // Validate direction across three successively extended time windows.
      out.ready=gravity_stable_windows_>=3;out.gravity_estimated=out.ready;
      return out;
    }
    auto fitted=fit(observations);out.fit_error=fitted.error;
    if(!fitted.valid)return out;
    out.ready=true;out.position=translation;out.velocity=fitted.initial+dv;
    out.rotation=fused_rotation;out.gravity=gravity;out.velocity_cov=fitted.covariance;
    return out;
  }
private:
  PoseMatrix last_pose_covariance=PoseMatrix::Zero();
  M3D last_integration_correction_=M3D::Identity();
  V3D bg_=V3D::Zero(),ba_=V3D::Zero(),previous_omega_=V3D::Zero();
  double scale_=1.;
  M3D gyro_rotation_=M3D::Identity();
  struct InertialSample { double time; M3D rotation; V3D acceleration; };
  std::vector<InertialSample> pending_imus;
  InertialSample last_sample{0,M3D::Identity(),V3D::Zero()};
  M3D last_rotation=M3D::Identity(),last_imu_rotation=M3D::Identity();
  void integrateAcceptedMotion(const M3D &correction,const V3D &gravity,double end) {
    const Eigen::Quaterniond first(last_integration_correction_);
    const Eigen::Quaterniond last(correction);
    double cursor=last_scan_time;
    auto advance=[&](const InertialSample &sample,double time) {
      if(time<=cursor)return;
      const double alpha=std::max(0.,std::min(1.,(time-last_scan_time)/(end-last_scan_time)));
      const V3D acc=first.slerp(alpha,last)*(sample.rotation*sample.acceleration)+gravity;
      const double dt=time-cursor;const V3D average=.5*(previous_acc+acc);
      dp+=dv*dt+.5*average*dt*dt;dv+=average*dt;
      previous_acc=acc;cursor=time;
    };
    for(const auto &sample:pending_imus) {advance(sample,sample.time);last_sample=sample;}
    if(std::abs(cursor-end)>1e-7)throw std::runtime_error("motion initialization requires an IMU endpoint bracket");
    pending_imus.clear();
  }
  bool timeout_fallback_;
  double fallback_velocity_sigma_;
  bool estimate_gravity_;
  int gravity_stable_windows_=0;
  V3D last_fit_gravity_=V3D::Zero();
  bool setTarget(const PointCloudXYZI::Ptr &cloud) {
    target=cloud;
    grid.clear();
    buildGrid();
    normals.assign(target->size(),V3D::Zero());
    centers.resize(target->size());
    size_t usable=0;
    for(size_t i=0;i<target->size();++i) {
      std::vector<int> ids(12);std::vector<float> dist(12);
      if(nearest((*target)[i],12,ids,dist)!=12 || dist.back()>2.25)continue;
      V3D mean=V3D::Zero();for(int id:ids)mean+=(*target)[id].getVector3fMap().cast<double>();mean/=12.;
      M3D cov=M3D::Zero();for(int id:ids){V3D d=(*target)[id].getVector3fMap().cast<double>()-mean;cov+=d*d.transpose();}
      Eigen::SelfAdjointEigenSolver<M3D> es(cov/12.);
      if(es.eigenvalues()[1]<.005 || es.eigenvalues()[0]>.01 || es.eigenvalues()[0]>.15*es.eigenvalues()[1])continue;
      normals[i]=es.eigenvectors().col(0);centers[i]=mean;++usable;
    }
    return usable>=100;
  }
  void resetWindow() {
    start=previous_time=last_scan_time=-1;target.reset();grid.clear();normals.clear();centers.clear();observations.clear();
    previous_acc.setZero();dp.setZero();dv.setZero();velocity_guess.setZero();last_translation.setZero();
    pending_imus.clear();last_rotation.setIdentity();last_imu_rotation.setIdentity();
    last_pose_covariance.setZero();last_integration_correction_.setIdentity();
    gravity_stable_windows_=0;last_fit_gravity_.setZero();
  }
  double first_attempt=-1,start=-1,previous_time=-1,last_scan_time=-1;
  Eigen::Quaterniond anchor=Eigen::Quaterniond::Identity();
  V3D previous_acc=V3D::Zero(),dp=V3D::Zero(),dv=V3D::Zero(),velocity_guess=V3D::Zero(),last_translation=V3D::Zero();
  PointCloudXYZI::Ptr target;
  struct Hash {size_t operator()(const std::array<int,3>& a)const {return size_t(a[0])*73856093u ^ size_t(a[1])*19349663u ^ size_t(a[2])*83492791u;}};
  std::unordered_map<std::array<int,3>,std::vector<int>,Hash> grid;
  static std::array<int,3> gridKey(const PointType &p){return {{int(std::floor(p.x/.8)),int(std::floor(p.y/.8)),int(std::floor(p.z/.8))}};}
  void buildGrid(){for(size_t i=0;i<target->size();++i)grid[gridKey((*target)[i])].push_back(i);}
  int nearest(const PointType &p,int k,std::vector<int>& ids,std::vector<float>& distances)const {
    auto key=gridKey(p);int radius=k==1?1:2;std::vector<std::pair<float,int>> found;
    for(int x=-radius;x<=radius;++x)for(int y=-radius;y<=radius;++y)for(int z=-radius;z<=radius;++z){
      auto it=grid.find({{key[0]+x,key[1]+y,key[2]+z}});if(it==grid.end())continue;
      for(int id:it->second){const auto &q=(*target)[id];float d=(p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y)+(p.z-q.z)*(p.z-q.z);if(d<=(k==1?.64f:2.25f))found.emplace_back(d,id);}
    }
    int n=std::min(k,int(found.size()));std::partial_sort(found.begin(),found.begin()+n,found.end());
    for(int i=0;i<n;++i){distances[i]=found[i].first;ids[i]=found[i].second;}return n;
  }
  std::vector<V3D> normals,centers;
  std::vector<Observation> observations;
};
