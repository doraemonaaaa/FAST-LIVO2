#include <iostream>
#include "initialization/batch_imu_initializer.hpp"
#include "initialization/imu_time_offset/lidar_rotation.hpp"
#include <stdexcept>

namespace initialization {
BatchLidarTrajectory BuildBatchLidarTrajectory(const std::vector<LivoScan>& scans,
    const Eigen::Matrix4d& imu_lidar, int workers) {
    using namespace imu_time_offset;
    BatchLidarTrajectory result;
    std::vector<RotationCloud> raw;
    double first=-1;
    for(const auto& scan:scans) {
        if(scan.end-scan.begin<.09) continue;
        if(first<0) first=scan.begin;
        if(scan.begin-first>30.) break;
        raw.push_back(FullRotationCloud(RotationClouds(scan,imu_lidar)));
        result.observation_end=scan.end;
    }
    result.scans=raw.size();
    if(raw.size()<20) throw std::runtime_error("batch IMU initialization needs at least 20 full LiDAR scans");
    std::vector<RotationFit> fits(raw.size());auto corrected=raw;
    for(int pass=0;pass<3;++pass) {
        auto next=RegisterRotationPairs(corrected,fits,workers);
        size_t valid=0;
        for(size_t i=1;i<raw.size();++i) valid+=next[i].valid;
        fits=std::move(next);
        std::cout<<"[batch initialization] LiDAR-only pass "<<pass<<": "<<valid<<'/'<<raw.size()-1<<" valid pairs";
        if(pass==2) break;
        // Fits are frozen; each task owns its corrected scan. The task barrier
        // finishes every deskew before the next registration pass starts.
        ForEachInitializationTask(raw.size()-2,workers,[&](size_t job) {
            const size_t i=job+1;
            const auto twist=EstimateTwist(fits[i],raw[i].time-raw[i-1].time,
                                          fits[i+1],raw[i+1].time-raw[i].time);
            corrected[i]=twist.valid ? DeskewRotationCloud(raw[i],twist) : raw[i];
        });
    }
    // Failed registrations split the trajectory. A predicted displacement is
    // never silently promoted to a LiDAR observation.
    std::vector<BatchLidarPose> segment{{raw[0].time,V3D::Zero(),M3D::Identity()}},best;
    auto retain=[&]() {
        if(best.empty() || segment.back().time-segment.front().time>best.back().time-best.front().time)
            best=segment;
    };
    for(size_t i=1;i<raw.size();++i) {
        if(!fits[i].valid) {
            ++result.rejected_pairs;retain();
            segment={{raw[i].time,V3D::Zero(),M3D::Identity()}};
            continue;
        }
        const auto& previous=segment.back();
        segment.push_back({raw[i].time,previous.position+previous.rotation*fits[i].translation,
                           previous.rotation*fits[i].rotation});
    }
    retain();
    if(best.back().time-best.front().time<3.)
        throw std::runtime_error("batch IMU initialization found no 3-second continuous LiDAR trajectory");
    // Identical to the reviewed experiment: every fifth scan (~0.5 s), first
    // 10 seconds of the longest accepted segment, actual sensor timestamps.
    for(size_t i=0;i<best.size();i+=5) {
        if(best[i].time-best.front().time>10.) break;
        result.poses.push_back(best[i]);
    }
    if(result.poses.size()<6)
        throw std::runtime_error("insufficient batch calibration keyframes");
    return result;
}
}
