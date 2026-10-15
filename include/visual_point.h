/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#ifndef LIVO_POINT_H_
#define LIVO_POINT_H_

#include <boost/noncopyable.hpp>
#include "common_lib.h"
#include "frame.h"

class Feature;

/// A visual map point on the surface of the scene.
class VisualPoint : boost::noncopyable
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Vector3d pos_;                //!< 3d pos of the point in the world coordinate frame.
  Vector3d normal_;             //!< Surface normal at point.
  Matrix3d normal_information_; //!< Inverse covariance matrix of normal estimation.
  Vector3d previous_normal_;    //!< Last updated normal vector.
  struct CameraObservation
  {
    list<Feature *> features;
    Feature *reference = nullptr;
    bool has_reference = false;
    bool reference_converged = false;
  };
  vector<CameraObservation> camera_observations;
  list<Feature *> &observations(size_t camera) { return camera_observations.at(camera).features; }
  const list<Feature *> &observations(size_t camera) const { return camera_observations.at(camera).features; }
  Feature *&referencePatch(size_t camera) { return camera_observations.at(camera).reference; }
  bool &hasReference(size_t camera) { return camera_observations.at(camera).has_reference; }
  bool &referenceConverged(size_t camera) { return camera_observations.at(camera).reference_converged; }
  Eigen::Matrix3d covariance_;  //!< Covariance of the point.
  bool is_normal_initialized_;  //!< True if the normal is initialized.

  VisualPoint(const Vector3d &pos, size_t camera_count);
  ~VisualPoint();
  void findMinScoreFeature(const Vector3d &framepos, Feature *&ftr, size_t camera) const;
  void deleteNonRefPatchFeatures(size_t camera);
  void deleteFeatureRef(Feature *ftr);
  void addFrameRef(Feature *ftr);
  bool getCloseViewObs(const Vector3d &pos, Feature *&obs, const Vector2d &cur_px, size_t camera) const;
};

#endif // LIVO_POINT_H_
