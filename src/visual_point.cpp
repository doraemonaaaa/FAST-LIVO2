/* 
This file is part of FAST-LIVO2: Fast, Direct LiDAR-Inertial-Visual Odometry.

Developer: Chunran Zheng <zhengcr@connect.hku.hk>

For commercial use, please contact me at <zhengcr@connect.hku.hk> or
Prof. Fu Zhang at <fuzhang@hku.hk>.

This file is subject to the terms and conditions outlined in the 'LICENSE' file,
which is included as part of this source code package.
*/

#include "visual_point.h"
#include "feature.h"
#include <stdexcept>
#include <vikit/math_utils.h>

VisualPoint::VisualPoint(const Vector3d &pos, size_t camera_count)
    : pos_(pos), previous_normal_(Vector3d::Zero()), normal_(Vector3d::Zero()),
      is_normal_initialized_(false), camera_observations(camera_count)
{
}

VisualPoint::~VisualPoint()
{
  for (auto &camera : camera_observations)
    for (auto *feature : camera.features) delete feature;
}

void VisualPoint::addFrameRef(Feature *ftr)
{
  observations(ftr->camera_index).push_front(ftr);
}

void VisualPoint::deleteFeatureRef(Feature *ftr)
{
  auto &camera = camera_observations.at(ftr->camera_index);
  if (camera.reference == ftr) { camera.reference = nullptr; camera.has_reference = false; }
  auto &features = camera.features;
  for (auto it = features.begin(); it != features.end(); ++it)
    if (*it == ftr) { delete *it; features.erase(it); return; }
}

bool VisualPoint::getCloseViewObs(const Vector3d &framepos, Feature *&ftr, const Vector2d &cur_px, size_t camera) const
{
  // TODO: get frame with same point of view AND same pyramid level!
  if (observations(camera).size() <= 0) return false;

  Vector3d obs_dir(framepos - pos_);
  obs_dir.normalize();
  auto min_it = observations(camera).begin();
  double min_cos_angle = 0;
  for (auto it = observations(camera).begin(), ite = observations(camera).end(); it != ite; ++it)
  {
    Vector3d dir((*it)->T_f_w_.inverse().translation() - pos_);
    dir.normalize();
    double cos_angle = obs_dir.dot(dir);
    if (cos_angle > min_cos_angle)
    {
      min_cos_angle = cos_angle;
      min_it = it;
    }
  }
  ftr = *min_it;

  // Vector2d ftr_px = ftr->px_;
  // double pixel_dist = (cur_px-ftr_px).norm();

  // if(pixel_dist > 200)
  // {
  //   ROS_ERROR("The pixel dist exceeds 200.");
  //   return false;
  // }

  if (min_cos_angle < 0.5) // assume that observations larger than 60° are useless 0.5
  {
    // ROS_ERROR("The obseved angle is larger than 60°.");
    return false;
  }

  return true;
}

void VisualPoint::findMinScoreFeature(const Vector3d &framepos, Feature *&ftr, size_t camera) const
{
  auto min_it = observations(camera).begin();
  float min_score = std::numeric_limits<float>::max();

  for (auto it = observations(camera).begin(), ite = observations(camera).end(); it != ite; ++it)
  {
    if ((*it)->score_ < min_score)
    {
      min_score = (*it)->score_;
      min_it = it;
    }
  }
  ftr = *min_it;
}

void VisualPoint::deleteNonRefPatchFeatures(size_t camera)
{
  for (auto it = observations(camera).begin(); it != observations(camera).end();)
  {
    if (*it != camera_observations.at(camera).reference)
    {
      delete *it;
      it = observations(camera).erase(it);
    }
    else
    {
      ++it;
    }
  }
}