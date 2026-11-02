#pragma once
#include <Eigen/Core>

// A pending scan has not entered the map yet. Keep its deskewed coordinates
// relative to the current IMU frame when a visual update corrects that frame.
// Applying T_after * inverse(T_before) to every pending world point avoids
// mixing pieces expressed using different posterior poses in the next LIO.
template<class Points>
void correctPendingCloud(Points &points, const Eigen::Matrix3d &before_rotation,
    const Eigen::Vector3d &before_position, const Eigen::Matrix3d &after_rotation,
    const Eigen::Vector3d &after_position)
{
  const Eigen::Matrix3d rotation = after_rotation * before_rotation.transpose();
  for (auto &sample : points)
    sample.world = rotation * (sample.world - before_position) + after_position;
}
