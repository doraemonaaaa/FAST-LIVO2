#pragma once
#include "common_lib.h"
#include <map>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <Eigen/Eigenvalues>
#include <deque>

// Short-window, attitude-aided LiDAR translation bootstrap. No map/EKF feedback.
// Points and inertial increments are expressed in the first scan's IMU frame.
class MotionInitializer {
public:
  explicit MotionInitializer(bool timeout_fallback=false, double velocity_sigma=1.0)
      : timeout_fallback_(timeout_fallback), fallback_velocity_sigma_(std::max(1.0, velocity_sigma)) {}
  struct Result {
    bool ready=false, timeout_fallback=false;
    V3D velocity=V3D::Zero(), position=V3D::Zero(), gravity=V3D::Zero();
    M3D rotation=M3D::Identity(), velocity_cov=M3D::Identity();
    double elapsed=0, match_ratio=0, residual=0, eigen_ratio=0, fit_error=0;
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
  Result update(const PointCloudXYZI &scan, double begin, double end,
                const std::deque<sensor_msgs::Imu::ConstPtr> &imus,
                const M3D &R_LI, const V3D &t_LI) {
    Result out;
    if(imus.empty() || scan.empty()) return out;
    const auto &m=imus.back(); Eigen::Quaterniond q(m->orientation.w,m->orientation.x,m->orientation.y,m->orientation.z);
    if(m->orientation_covariance[0]<0 || !q.coeffs().allFinite() || std::abs(q.norm()-1)>.01)
      throw std::runtime_error("motion initialization: invalid IMU attitude");
    q.normalize();
    if(first_attempt<0)first_attempt=end;
    if(start<0) { start=end; anchor=q; previous_time=end; previous_acc=V3D::Zero(); }
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
    M3D R=(anchor.conjugate()*q).toRotationMatrix();
    V3D gravity=anchor.conjugate()*V3D(0,0,-9.81);
    // Integrate IMU acceleration at its own timestamps; do not use scan-rate acceleration.
    for(const auto &im:imus) {
      double time=im->header.stamp.toSec(); if(time<=previous_time) continue;
      double dt=time-previous_time;
      if(dt>.1) {
        ROS_WARN("[motion init] IMU gap %.6f s: discard startup window and restart",dt);
        resetWindow();
        return update(scan,begin,end,imus,R_LI,t_LI);
      }
      Eigen::Quaterniond qi(im->orientation.w,im->orientation.x,im->orientation.y,im->orientation.z);
      if(!qi.coeffs().allFinite() || std::abs(qi.norm()-1)>.01) throw std::runtime_error("motion initialization: invalid attitude in window");
      V3D acc=(anchor.conjugate()*qi.normalized())*V3D(im->linear_acceleration.x,im->linear_acceleration.y,im->linear_acceleration.z)+gravity;
      V3D av=.5*(previous_acc+acc); dp+=dv*dt+.5*av*dt*dt; dv+=av*dt; previous_acc=acc;previous_time=time;
    }
    // End-of-scan attitude, rotational deskew, and the latest independent velocity estimate.
    V3D omega(m->angular_velocity.x,m->angular_velocity.y,m->angular_velocity.z);
    PointCloudXYZI::Ptr cloud(new PointCloudXYZI());
    for(const auto &pt:scan) {
      V3D p(pt.x,pt.y,pt.z); if(!p.allFinite() || p.norm()<1 || p.norm()>45) continue;
      double dt=end-(begin+double(pt.curvature)*.001);
      if(dt<-.01 || dt>.3) continue;
      V3D w=-omega*dt; double angle=w.norm();M3D dR=M3D::Identity();if(angle>1e-10)dR=Eigen::AngleAxisd(angle,w/angle).toRotationMatrix();
      p=R*dR*(R_LI*p+t_LI)-velocity_guess*dt;
      PointType v=pt;v.x=p.x();v.y=p.y();v.z=p.z();cloud->push_back(v);
    }
    // Deterministic voxel centroids; avoid binary PCL/Eigen alignment dependencies.
    std::map<std::array<int,3>,std::pair<V3D,int>> cells;
    for(const auto &p:*cloud) {
      std::array<int,3> key{{int(std::floor(p.x/.25)),int(std::floor(p.y/.25)),int(std::floor(p.z/.25))}};
      auto it=cells.find(key);
      if(it==cells.end())it=cells.emplace(key,std::make_pair(V3D::Zero().eval(),0)).first;
      it->second.first+=V3D(p.x,p.y,p.z);++it->second.second;
    }
    PointCloudXYZI::Ptr reduced(new PointCloudXYZI());
    for(const auto &cell:cells) {V3D q=cell.second.first/cell.second.second;PointType p;p.x=q.x();p.y=q.y();p.z=q.z();reduced->push_back(p);}
    if(reduced->size()<100) return out;
    if(!target) { target=reduced; buildGrid();buildNormals();return out; }
    V3D translation=last_translation+velocity_guess*(end-last_scan_time);
    if(last_scan_time<0)translation.setZero();
    M3D H=M3D::Zero();double rms=0;int count=0;
    for(int iteration=0;iteration<20;++iteration) {
      H.setZero();V3D gradient=V3D::Zero();double squared=0;count=0;
      for(const auto &pt:*reduced) {
        V3D p=V3D(pt.x,pt.y,pt.z)+translation;PointType search;search.x=p.x();search.y=p.y();search.z=p.z();
        std::vector<int> ids(1);std::vector<float> distances(1);
        if(nearest(search,1,ids,distances)!=1 || distances[0]>.64 || normals[ids[0]].norm()<.5)continue;
        const V3D &n=normals[ids[0]];double r=n.dot(p-centers[ids[0]]);if(std::abs(r)>.3)continue;
        double weight=std::min(1.,.05/std::max(1e-9,std::abs(r)));
        H+=weight*n*n.transpose();gradient+=weight*n*r;squared+=r*r;++count;
      }
      if(count<100)return out;
      Eigen::SelfAdjointEigenSolver<M3D> es(H);out.eigen_ratio=es.eigenvalues().minCoeff()/es.eigenvalues().maxCoeff();
      if(out.eigen_ratio<.005)return out;
      V3D step=-H.ldlt().solve(gradient);if(!step.allFinite() || step.norm()>.75)return out;
      translation+=step;rms=std::sqrt(squared/count);if(step.norm()<1e-4)break;
    }
    out.match_ratio=double(count)/reduced->size();out.residual=rms;
    if(out.match_ratio<.20 || rms>.10)return out;
    // Use isotropic conservative per-scan registration noise with observed anisotropy.
    M3D info=H/double(count)/(.03*.03);
    observations.push_back({out.elapsed,translation,dp,dv,info});
    last_translation=translation;last_scan_time=end;
    if(out.elapsed>.15)velocity_guess=(translation-dp)/out.elapsed+dv;
    if(out.elapsed<.8)return out;
    auto fitted=fit(observations);out.fit_error=fitted.error;
    if(!fitted.valid)return out;
    out.ready=true;out.position=translation;out.velocity=fitted.initial+dv;
    out.rotation=R;out.gravity=gravity;out.velocity_cov=fitted.covariance;
    return out;
  }
private:
  bool timeout_fallback_;
  double fallback_velocity_sigma_;
  void buildNormals() {
    normals.resize(target->size(),V3D::Zero());centers.resize(target->size(),V3D::Zero());
    for(size_t i=0;i<target->size();++i) {
      std::vector<int> ids(12);std::vector<float> dist(12);
      if(nearest((*target)[i],12,ids,dist)!=12 || dist.back()>2.25)continue;
      V3D mean=V3D::Zero();for(int id:ids)mean+=(*target)[id].getVector3fMap().cast<double>();mean/=12.;
      M3D cov=M3D::Zero();for(int id:ids){V3D d=(*target)[id].getVector3fMap().cast<double>()-mean;cov+=d*d.transpose();}
      Eigen::SelfAdjointEigenSolver<M3D> es(cov/12.);
      if(es.eigenvalues()[1]<.005 || es.eigenvalues()[0]>.01 || es.eigenvalues()[0]>.15*es.eigenvalues()[1])continue;
      normals[i]=es.eigenvectors().col(0);centers[i]=mean;
    }
  }
  void resetWindow() {
    start=previous_time=last_scan_time=-1;target.reset();grid.clear();normals.clear();centers.clear();observations.clear();
    previous_acc.setZero();dp.setZero();dv.setZero();velocity_guess.setZero();last_translation.setZero();
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
