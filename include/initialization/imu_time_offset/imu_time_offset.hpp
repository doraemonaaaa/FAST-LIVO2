#pragma once

#include "initialization/imu_time_offset/estimator.hpp"
#include "initialization/input.h"
#include <filesystem>
#include <json/json.h>

namespace initialization {

struct OffsetConfig {
    double imu_time_offset_s=0, initialization_acceleration_norm=9.81;
};
// One immutable configured-clock state per temporary LIVO update. All offset
// candidates use the same interpolated state at each query's reference epoch.
struct OffsetInitializationState {
    double time = 0;
    Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity();
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero(), gravity = Eigen::Vector3d::Zero();
    Eigen::Vector3d gyro_bias = Eigen::Vector3d::Zero(), acceleration_bias = Eigen::Vector3d::Zero();
};

livo::ImuTimeOffsetResult EstimateImuTimeOffset(
    const std::vector<LivoScan>& scans,
    const std::vector<ImuSample::ConstPtr>& imu,
    const OffsetConfig& param,
    const Eigen::Matrix4d& imu_lidar,
    const std::filesystem::path& diagnostics,
    double window_begin, double window_end,
    const std::vector<OffsetInitializationState>& states);

Json::Value ImuTimeOffsetJson(const livo::ImuTimeOffsetResult& result,
                             double configured_offset_s);
}  // namespace initialization
