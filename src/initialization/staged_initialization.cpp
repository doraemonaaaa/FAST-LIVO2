#include "initialization/staged_initialization.h"
#include "LIVMapper.h"
#include <iomanip>
#include <stdexcept>

void StagedInitialization::configure(ros::NodeHandle& nh) {
  nh.param<bool>("imu/staged_initialization",enabled_,false);
  nh.param<int>("imu/initialization_threads",workers_,initialization::kDefaultInitializationThreads);
  if(workers_<1)throw std::runtime_error("imu/initialization_threads must be positive");
  nh.param<std::string>("imu/staged_phase",phase_,"refine");
  nh.param<std::string>("imu/staged_seed_output",output_,"");
  if(enabled_ && (phase_!="window" && phase_!="refine"))
    throw std::runtime_error("Stage one must use the independent batch initializer");
  if(enabled_ && phase_=="window") {
    if(output_.empty())throw std::runtime_error("Missing stage-two clock output");
    states_.open(output_+".states.csv");
    if(!states_)throw std::runtime_error("Cannot create staged state audit");
    states_<<"time,bgx,bgy,bgz,bax,bay,baz,vx,vy,vz,gx,gy,gz,bg_var_x,bg_var_y,bg_var_z,ba_var_x,ba_var_y,ba_var_z\n";
  }
}
void StagedInitialization::observeScan(const PointCloudXYZI& points,double begin) {
  if(enabled_ && phase_=="window" && !points.empty())
    scans_.push_back(initialization::CalibrationScan(points,begin));
}
void StagedInitialization::observe(const LidarMeasureGroup& group) {
  if(!enabled_ || phase_!="window")return;
  for(const auto& sample:group.measures.back().imu) {
    auto imu=initialization::CalibrationImu(*sample);
    if(imus_.empty() || imu->timestamp>imus_.back()->timestamp)imus_.push_back(imu);
  }
}
bool StagedInitialization::advance(LIVMapper& mapper) {
  if(!enabled_ || phase_=="refine")return true;
  const auto& group=mapper.LidarMeasures;const auto& measure=group.measures.back();
  const double time=group.lio_vio_flg==LIO?measure.lio_time:measure.vio_time;
  if(time<=last_)return false;
  if(ready_<0) {if(group.lio_vio_flg!=LIO)return false;ready_=time;}
  last_=time;
  const auto& state=mapper._state;
  snapshots_.push_back({time,state.rot_end,state.vel_end,state.gravity,state.bias_g,state.bias_a});
  states_<<std::setprecision(17)<<time;
  for(const auto* v:{&state.bias_g,&state.bias_a,&state.vel_end,&state.gravity})
    for(int j=0;j<3;++j)states_<<','<<(*v)[j];
  for(int j=10;j<16;++j)states_<<','<<state.cov(j,j);
  states_<<'\n'<<std::flush;
  ROS_INFO("[staged init] post-motion clock window %.3f / 3.000 s",time-ready_);
  if(time-ready_<3.)return false;
  const double delta=chooseOffset(mapper);
  Json::Value seed;seed["imu_offset_delta"]=delta;
  seed["configured_imu_time_offset_s"]=mapper.imu_time_offset;
  seed["effective_imu_time_offset_s"]=mapper.imu_time_offset+delta;
  seed["bias_handoff"]="retain batch means; recalibrate if the clock changes";
  for(int i=0;i<3;++i) {
    seed["temporary_posterior_bias_g"].append(state.bias_g[i]);
    seed["temporary_posterior_bias_a"].append(state.bias_a[i]);
  }
  std::ofstream out(output_);out<<seed<<'\n';out.flush();
  if(!out)throw std::runtime_error("Cannot save stage-two clock decision");
  return true;
}
double StagedInitialization::chooseOffset(LIVMapper& mapper) {
  // Use a small epoch for the same preintegration arithmetic as the reference.
  const double epoch=std::floor(ready_)-1.;auto scans=scans_;auto snapshots=snapshots_;
  for(auto& scan:scans) {scan.begin-=epoch;scan.end-=epoch;for(auto& p:scan.points)p.time-=epoch;}
  for(auto& state:snapshots)state.time-=epoch;
  std::vector<initialization::ImuSample::ConstPtr> imu;
  for(const auto& sample:imus_) {auto p=std::make_shared<initialization::ImuSample>(*sample);p->timestamp-=epoch;imu.push_back(p);}
  Eigen::Matrix4d extrinsic=Eigen::Matrix4d::Identity();extrinsic.topLeftCorner<3,3>()=mapper.extR;extrinsic.topRightCorner<3,1>()=mapper.extT;
  initialization::OffsetConfig config;config.imu_time_offset_s=mapper.imu_time_offset;
  config.workers=workers_;
  config.initialization_acceleration_norm=mapper.p_imu->IMU_mean_acc_norm;
  config.progress=[&mapper]{mapper.startup_replay.reportProgress();};
  const auto directory=std::filesystem::path(output_).parent_path();
  auto result=initialization::EstimateImuTimeOffset(scans,imu,config,extrinsic,directory,
    ready_-epoch,last_-epoch,snapshots);
  std::ofstream(directory/"imu_time_offset.json")<<initialization::ImuTimeOffsetJson(result,config.imu_time_offset_s)<<'\n';
  ROS_WARN("[staged init] clock applied=%d stamp shift %.6f s: %s",int(result.applied),result.stamp_shift_s,result.reason);
  return -result.stamp_shift_s; // Callback subtracts imu_time_offset.
}
