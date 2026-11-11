#pragma once
#include "initialization/input.h"
#include <filesystem>
#include <json/json.h>

namespace initialization {


// LiDAR-only poses in the first accepted scan's IMU frame. Registration
// uncertainty is represented by relative-pose factors in the batch solve.
struct BatchLidarPose {
    double time=0;
    V3D position=V3D::Zero();
    M3D rotation=M3D::Identity();
};
struct BatchLidarTrajectory {
    std::vector<BatchLidarPose> poses;
    double observation_end=0;
    size_t scans=0, rejected_pairs=0;
};
struct BatchImuInitialization {
    V3D gyro_bias=V3D::Zero(), acceleration_bias=V3D::Zero();
    V3D gravity=V3D(0,0,-9.81), initial_velocity=V3D::Zero();
    double calibration_begin=0, calibration_end=0;
    Json::Value report;
};

BatchLidarTrajectory BuildBatchLidarTrajectory(const std::vector<LivoScan>& scans,
    const Eigen::Matrix4d& imu_lidar);

// Bias is approximately constant within this short calibration window. Both
// biases and gravity remain estimated by the normal ESIKF after handoff.
BatchImuInitialization SolveBatchImuInitialization(
    const std::vector<BatchLidarPose>& poses,
    const std::vector<ImuSample::ConstPtr>& imu,
    const std::filesystem::path& diagnostics={});
}
