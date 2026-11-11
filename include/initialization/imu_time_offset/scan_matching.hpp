#pragma once

#include "initialization/imu_time_offset/lidar_rotation.hpp"
#include "initialization/imu_time_offset/estimator.hpp"
#include <memory>

namespace initialization { namespace imu_time_offset {

// Initial state is identical for every candidate of a query. Velocity comes
// from independent LiDAR registration; gravity and biases stay frozen.
struct DeskewState {
    V3D velocity = V3D::Zero(), gravity = V3D::Zero();
    V3D gyro_bias = V3D::Zero(), acceleration_bias = V3D::Zero();
    double acceleration_scale = 1.;
};

bool HasImuCoverage(const RotationCloud& input,
                    const std::vector<ImuSample::ConstPtr>& imu,
                    double padding);
RotationCloud DeskewWithImu(const RotationCloud& input,
                            const std::vector<ImuSample::ConstPtr>& imu,
                            const DeskewState& state, double stamp_shift_s);
DeskewState QueryDeskewState(const RotationCloud& input,
                             const std::vector<ImuSample::ConstPtr>& imu,
                             const LidarTwist& twist, double acceleration_scale);

// Owns immutable target points, normals and spatial index. Registration starts
// from identity for every offset and rescans associations at the final pose.
class FrozenPlaneMap {
 public:
    explicit FrozenPlaneMap(const RotationCloud& target);
    ~FrozenPlaneMap();
    FrozenPlaneMap(FrozenPlaneMap&&) noexcept;
    FrozenPlaneMap& operator=(FrozenPlaneMap&&) noexcept;
    livo::OffsetMatchScore Evaluate(const RotationCloud& source) const;
 private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} }  // namespace initialization::imu_time_offset
