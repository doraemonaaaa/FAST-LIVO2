/* Parameterized visual observations and one joint iterated EKF update. */
#include "vio.h"
#include "state_numerics.h"
#include <unordered_set>
#include <stdexcept>
#include <iomanip>

VIOManager::~VIOManager()
{
  cameras.clear();
  for (auto &entry : visual_map) delete entry.second;
}

void VIOManager::clearStartupMap()
{
  for (auto &camera : cameras) camera->clearStartupMap();
  for (auto &entry : visual_map) delete entry.second;
  visual_map.clear();
}

void VIOManager::processFrame(const vector<MeasureGroup::Image> &images,
    vector<pointWithVar> &pg, const unordered_map<VOXEL_LOCATION, VoxelOctoTree *> &map)
{
  for (auto &camera : cameras) camera->active = false;
  for (const auto &image : images)
  {
    auto &camera = *cameras.at(image.camera_index);
    if (camera.active) throw std::runtime_error("Duplicate camera in image bundle");
    camera.active = true;
    camera.prepareFrame(image, pg, map);
  }
  // One geometric point contributes through one camera per bundle. This avoids
  // counting its shared depth/normal uncertainty as independent information.
  // Selection is deterministic, in configured camera priority order.
  std::unordered_set<VisualPoint *> claimed;
  for (auto &camera : cameras)
    if (camera->active)
      for (int i = 0; i < camera->total_points; ++i)
        if (!claimed.insert(camera->visual_submap->voxel_points[i]).second)
          camera->noise_weights[i] = 0;
  static const bool skip_visual = [] {
    const char *v = std::getenv("FLIVO_DIAGNOSTIC_SKIP_VISUAL");
    return v && string(v) == "1";
  }();
  if (!skip_visual)
  {
    computeJacobianAndUpdateEKF();
    // All views used the same pre-update map. Commit patches only afterwards.
    for (auto &camera : cameras)
      if (camera->active) camera->finishFrame(pg, map);
  }
  if (audit)
    for (const auto &camera : cameras)
      if (camera->active)
      {
        size_t used = 0;
        for (double w : camera->noise_weights) if (w != 0) ++used;
        audit << std::setprecision(17) << camera->image_timestamp << ',' << camera->id << ','
              << camera->total_points << ',' << used << ',' << state->exposure(camera->camera_index)
              << ',' << state->cov.trace() << '\n';
      }
}

void VIOManager::computeJacobianAndUpdateEKF()
{
  const int dimension = state->dimension();
  const int columns = 6 + cameras.size();
  vector<int> indices = {0, 1, 2, 3, 4, 5};
  for (size_t k = 0; k < cameras.size(); ++k) indices.push_back(StatesGroup::exposureIndex(k));
  int rows = 0, levels = 0, iterations = 0;
  for (auto &camera : cameras)
    if (camera->active)
    {
      rows += camera->total_points * camera->patch_size_total;
      levels = std::max(levels, camera->patch_pyrimid_level);
      iterations = std::max(iterations, camera->max_iterations);
    }
  if (!rows) return;
  G = MatrixXd::Zero(dimension, dimension);
  const double covariance_scale = anchor().img_point_cov;
  for (int level = levels - 1; level >= 0; --level)
  {
    G.setZero();
    StatesGroup old_state = *state;
    float last_error = std::numeric_limits<float>::max();
    for (auto &camera : cameras) camera->has_ref_patch_cache = false;
    for (int iteration = 0; iteration < iterations; ++iteration)
    {
      MatrixXd H = MatrixXd::Zero(rows, columns);
      VectorXd z = VectorXd::Zero(rows);
      float error = 0;
      int count = 0, offset = 0;
      for (auto &camera : cameras)
      {
        if (!camera->active) continue;
        const int n = camera->total_points * camera->patch_size_total;
        if (level >= camera->patch_pyrimid_level || !n) { offset += n; continue; }
        const int local_columns = camera->inverse_composition_en ? 6 : 7;
        MatrixXd local_H = MatrixXd::Zero(n, local_columns);
        VectorXd local_z = VectorXd::Zero(n);
        vector<float> errors(camera->total_points, 0);
        vector<int> counts(camera->total_points, 0);
        if (camera->inverse_composition_en)
          camera->buildInverseResiduals(camera->gray, level, local_H, local_z, errors, counts);
        else
          camera->buildResiduals(camera->gray, level, camera->noise_weights, local_H, local_z, errors, counts);
        const double weight = std::sqrt(covariance_scale / camera->img_point_cov);
        H.block(offset, 0, n, 6) = local_H.leftCols(6) * weight;
        if (local_columns == 7)
          H.block(offset, 6 + camera->camera_index, n, 1) = local_H.col(6) * weight;
        z.segment(offset, n) = local_z * weight;
        // Identical point/pixel accumulation order for N=1 and N>1.
        for (size_t i = 0; i < errors.size(); ++i)
        {
          error += float(errors[i] * weight * weight);
          count += counts[i];
        }
        offset += n;
      }
      if (!count) { *state = old_state; break; }
      error /= count;
      if (error > last_error) { *state = old_state; break; }
      old_state = *state;
      last_error = error;
      const MatrixXd information = H.transpose() * H;
      const VectorXd gradient = H.transpose() * z;
      MatrixXd full_information = MatrixXd::Zero(dimension, dimension);
      for (int i = 0; i < columns; ++i)
        for (int j = 0; j < columns; ++j)
          full_information(indices[i], indices[j]) = information(i, j);
      MatrixXd K1 = inverseStateMatrix(full_information + inverseStateMatrix(state->cov / covariance_scale));
      MatrixXd active_gain(dimension, columns);
      VectorXd prior_delta = *state_propagat - *state, active_delta(columns);
      for (int i = 0; i < columns; ++i)
      {
        active_gain.col(i) = K1.col(indices[i]);
        active_delta(i) = prior_delta(indices[i]);
      }
      MatrixXd gain_information = active_gain * information;
      gain_information.block<DIM_STATE,7>(0,0) = K1.block<DIM_STATE,7>(0,0)*information.block<7,7>(0,0);
      if (columns>7) gain_information.block<DIM_STATE,7>(0,0) += active_gain.block(0,7,DIM_STATE,columns-7)*information.block(7,0,columns-7,7);
      G.setZero();
      for (int i = 0; i < columns; ++i) G.col(indices[i]) = gain_information.col(i);
      VectorXd solution = -active_gain * gradient + prior_delta - gain_information * active_delta;
      solution.head<DIM_STATE>() = -K1.block<DIM_STATE,7>(0,0)*gradient.head<7>() + prior_delta.head<DIM_STATE>()
                                  - gain_information.block<DIM_STATE,7>(0,0)*active_delta.head<7>();
      if (columns>7) solution.head<DIM_STATE>() += -active_gain.block(0,7,DIM_STATE,columns-7)*gradient.tail(columns-7)
                                                -gain_information.block(0,7,DIM_STATE,columns-7)*active_delta.tail(columns-7);
      if (!solution.allFinite()) throw std::runtime_error("Non-finite joint visual update");
      *state += solution;
      if (solution.head<3>().norm() * 57.3f < 0.001f && solution.segment<3>(3).norm() * 100.f < 0.001f) break;
    }
  }
  const MD(DIM_STATE,DIM_STATE) covariance = state->cov.topLeftCorner<DIM_STATE,DIM_STATE>();
  const MD(DIM_STATE,DIM_STATE) base_gain = G.topLeftCorner<DIM_STATE,DIM_STATE>();
  MatrixXd posterior = state->cov-G*state->cov;
  posterior.topLeftCorner<DIM_STATE,DIM_STATE>() = covariance-base_gain*covariance;
  const int extra = dimension-DIM_STATE;
  if (extra) posterior.topLeftCorner<DIM_STATE,DIM_STATE>() -= G.topRightCorner(DIM_STATE,extra)*state->cov.bottomLeftCorner(extra,DIM_STATE);
  state->cov = posterior;
  for (auto &camera : cameras)
    if (camera->active) camera->updateFrameState(*state);
}
