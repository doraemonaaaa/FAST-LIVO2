#include "LIVMapper.h"
#include "initialization/gravity_norm_constraint.h"
#include "pending_cloud.h"

void LIVMapper::constrainGravityBeforePrediction(double target) {
  if(!gravity_norm_constraint || !gravity_est_en || startup_replay.pending() ||
     p_imu->imu_need_init || target<=p_imu->propagationTime())return;
  const M3D rotation=_state.rot_end;const V3D position=_state.pos_end;
  livo::ConstrainGravityNorm(_state,gravity_norm_sigma);
  // A correlated pose correction must not change pending scan body geometry.
  correctPendingCloud(inter_camera_points,rotation,position,_state.rot_end,_state.pos_end);
}
