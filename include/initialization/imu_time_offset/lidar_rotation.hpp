#pragma once

#include "initialization/input.h"
#include <array>

namespace initialization { namespace imu_time_offset {
constexpr int kPhases = 5;
constexpr double kSliceSeconds = .020;

struct RotationCloud {
    PointCloudXYZI::Ptr cloud{new PointCloudXYZI};
    std::vector<double> point_times;
    double time = 0;
};
using RotationPhases = std::array<RotationCloud,kPhases>;

struct RotationFit {
    bool valid = false;
    M3D rotation = M3D::Identity(), covariance = M3D::Identity();
    V3D translation = V3D::Zero();
    double ratio = 0, rms = 0;
    size_t matches = 0;
};

struct LidarTwist {
    bool valid = false;
    V3D angular_velocity = V3D::Zero(), velocity = V3D::Zero();
    V3D angular_acceleration = V3D::Zero(), acceleration = V3D::Zero();
};

RotationPhases RotationClouds(const LivoScan& scan, const Eigen::Matrix4d& transform);
RotationCloud FullRotationCloud(const RotationPhases& phases);
RotationFit RegisterRotation(const RotationCloud& target, const RotationCloud& source,
                             const RotationFit& prediction);
LidarTwist EstimateTwist(const RotationFit& backward, double dt_backward,
                         const RotationFit& forward, double dt_forward);
LidarTwist ShiftTwistReference(const LidarTwist& twist, double dt);
RotationCloud DeskewRotationCloud(const RotationCloud& cloud, const LidarTwist& twist);
} }  // namespace initialization::imu_time_offset
