#include "initialization/batch_imu_initializer.hpp"
#include "preprocess.h"
#include "imu_support.h"
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>

namespace {
template<class T>T Get(const YAML::Node& p,const char* key,T fallback) {
  return p[key]?p[key].as<T>():fallback;
}
}
// Offline, independent LiDAR/IMU calibration. No LIVMapper or ESIKF is created.
int main(int argc,char** argv) {
  try {
    if(argc!=5)throw std::runtime_error("Usage: fastlivo_batch_initialize PROFILE BAG OUTPUT_DIR BAG_START_SECONDS");
    ros::init(argc,argv,"batch_imu_initialization",ros::init_options::NoSigintHandler);
    const auto config=YAML::LoadFile(argv[1]);const auto imu_config=config["imu"];
    const auto directory=std::filesystem::path(argv[3]);std::filesystem::create_directories(directory);
    const int workers=Get<int>(imu_config,"initialization_threads",initialization::kDefaultInitializationThreads);
    if(workers<1)throw std::runtime_error("imu/initialization_threads must be positive");
    const auto input_start=std::chrono::steady_clock::now();
    const double requested_start=std::stod(argv[4]);
    if(!std::isfinite(requested_start) || requested_start<0)throw std::runtime_error("Invalid bag start");
    const double imu_offset=Get<double>(config["time_offset"],"imu_time_offset",0.);
    const double lidar_offset=Get<double>(config["time_offset"],"lidar_time_offset",0.);
    const double scale=9.81/Get<double>(imu_config,"initialization_acceleration_norm",9.81);
    if(!std::isfinite(scale) || scale<=0)throw std::runtime_error("Invalid acceleration normalization");
    const auto numerics=config["numerics"];
    if(Get<bool>(numerics,"fixed_eigen_cache",false)) {
      Eigen::setCpuCacheSizes(Get<int>(numerics,"eigen_l1_cache_bytes",48*1024),
        Get<int>(numerics,"eigen_l2_cache_bytes",2*1024*1024),Get<int>(numerics,"eigen_l3_cache_bytes",33*1024*1024));
    }
    Eigen::setNbThreads(1);
    const auto pre=config["preprocess"];Preprocess preprocess;
    preprocess.set(Get<bool>(pre,"feature_extract_enabled",false),Get<int>(pre,"lidar_type",AVIA),
      Get<double>(pre,"blind",.01),Get<int>(pre,"point_filter_num",3));
    preprocess.N_SCANS=Get<int>(pre,"scan_line",6);preprocess.blind_sqr=preprocess.blind*preprocess.blind;
    Eigen::Matrix4d extrinsic=Eigen::Matrix4d::Identity();
    const auto rotation=config["extrin_calib"]["extrinsic_R"].as<std::vector<double>>();
    const auto translation=config["extrin_calib"]["extrinsic_T"].as<std::vector<double>>();
    if(rotation.size()!=9 || translation.size()!=3)throw std::runtime_error("Invalid LiDAR-to-IMU extrinsic");
    for(int i=0;i<3;++i) {for(int j=0;j<3;++j)extrinsic(i,j)=rotation[i*3+j];extrinsic(i,3)=translation[i];}
    if(!extrinsic.allFinite())throw std::runtime_error("Non-finite extrinsic");
    const auto lidar_topic=config["common"]["lid_topic"].as<std::string>();
    const auto imu_topic=config["common"]["imu_topic"].as<std::string>();
    rosbag::Bag bag(argv[2],rosbag::bagmode::Read);rosbag::View all(bag);
    const double bag_begin=all.getBeginTime().toSec(),epoch=std::floor(bag_begin)-1.;
    const ros::Time start=ros::Time().fromSec(bag_begin+requested_start);
    rosbag::View view(bag,rosbag::TopicQuery(std::vector<std::string>{lidar_topic,imu_topic}),start);
    std::vector<initialization::LivoScan> scans;std::vector<double> records;
    std::vector<initialization::ImuSample::ConstPtr> imu;
    double first=-1,last_record=-1;
    for(const auto& message:view) {
      if(last_record>=0 && message.getTime().toSec()>last_record+31.)break;
      if(auto sample=message.instantiate<sensor_msgs::Imu>()) {
        auto corrected=*sample;corrected.header.stamp.fromSec(sample->header.stamp.toSec()-imu_offset-epoch);
        auto reading=initialization::CalibrationImu(corrected,scale);
        if(!imu.empty() && reading->timestamp<=imu.back()->timestamp)throw std::runtime_error("Non-monotonic IMU input");
        imu.push_back(reading);continue;
      }
      PointCloudXYZI::Ptr cloud(new PointCloudXYZI());double begin;
      if(auto sample=message.instantiate<sensor_msgs::PointCloud2>()) {
        begin=sample->header.stamp.toSec()+lidar_offset-epoch;preprocess.process(sample,cloud);
      } else if(auto sample=message.instantiate<livox_ros_driver::CustomMsg>()) {
        begin=sample->header.stamp.toSec()-epoch;preprocess.process(sample,cloud);
      } else continue;
      auto scan=initialization::CalibrationScan(*cloud,begin);
      if(scan.end-scan.begin<.09)continue;
      if(first<0) {first=scan.begin;last_record=message.getTime().toSec();}
      if(scan.begin-first<=30.) {scans.push_back(std::move(scan));records.push_back(message.getTime().toSec());}
    }
    if(imu.size()<2 || scans.empty())throw std::runtime_error("No usable startup LiDAR/IMU data");
    std::ofstream coverage(directory/"imu_coverage.csv");
    coverage<<"begin,end,points,status,reason\n"<<std::setprecision(17);
    std::vector<initialization::LivoScan> run,best;
    size_t discarded=0;double first_covered=-1;
    auto retain=[&]() {
      if(!run.empty() && (best.empty() || run.back().end-run.front().begin>best.back().end-best.front().begin))best=run;
    };
    for(const auto& scan:scans) {
      const auto status=imu_support::Check(imu,scan.begin,scan.end,[](const auto& s){return s->timestamp;});
      coverage<<scan.begin<<','<<scan.end<<','<<scan.points.size()<<','
        <<(status==imu_support::Coverage::Covered?"accepted":"dropped")<<','<<imu_support::Reason(status)<<'\n';
      if(status!=imu_support::Coverage::Covered) {++discarded;retain();run.clear();continue;}
      if(first_covered<0)first_covered=scan.begin;
      if(!run.empty() && imu_support::Check(imu,run.back().begin,scan.end,[](const auto& s){return s->timestamp;})!=imu_support::Coverage::Covered) {retain();run.clear();}
      run.push_back(scan);
    }
    retain();scans=std::move(best);
    if(scans.empty())throw std::runtime_error("No fully IMU-covered startup LiDAR scans");
    const double replay_start=requested_start;
    const auto trajectory_start=std::chrono::steady_clock::now();
    const auto trajectory=initialization::BuildBatchLidarTrajectory(scans,extrinsic,workers);
    const auto solve_start=std::chrono::steady_clock::now();
    auto result=initialization::SolveBatchImuInitialization(trajectory.poses,imu,directory);
    const auto solve_end=std::chrono::steady_clock::now();
    Json::Value timing;timing["workers"]=workers;
    timing["input_preparation_s"]=std::chrono::duration<double>(trajectory_start-input_start).count();
    timing["lidar_trajectory_s"]=std::chrono::duration<double>(solve_start-trajectory_start).count();
    timing["single_thread_bias_solve_s"]=std::chrono::duration<double>(solve_end-solve_start).count();
    std::ofstream(directory/"batch_initialization_timing.json")<<timing<<'\n';
    auto& report=result.report;report["epoch_s"]=epoch;report["scans"]=Json::UInt64(trajectory.scans);
    report["rejected_lidar_pairs"]=Json::UInt64(trajectory.rejected_pairs);
    report["observation_end_s"]=trajectory.observation_end;
    report["imu_time_offset_s"]=imu_offset;report["replay_bag_start_s"]=replay_start;
    report["uncovered_scans_discarded"]=Json::UInt64(discarded);
    report["first_covered_lidar_begin_s"]=first_covered;
    std::ofstream(directory/"batch_imu_initialization.json")<<report<<'\n';
    Json::Value seed;for(int i=0;i<3;++i) {seed["bias_g"].append(result.gyro_bias[i]);seed["bias_a"].append(result.acceleration_bias[i]);}
    seed["imu_offset_delta"]=0.;seed["replay_bag_start_s"]=replay_start;
    std::ofstream out(directory/"seed.yaml");out<<seed<<'\n';out.flush();
    if(!out)throw std::runtime_error("Cannot save batch bias means");
    std::cout<<"Batch bias ready; gyro "<<result.gyro_bias.transpose()<<" acceleration "<<result.acceleration_bias.transpose()<<'\n';
    return 0;
  } catch(const std::exception& e) {std::cerr<<"Batch initialization failed: "<<e.what()<<'\n';return 1;}
}
