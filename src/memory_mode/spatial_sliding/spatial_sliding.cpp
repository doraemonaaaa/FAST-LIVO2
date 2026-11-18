// Spatial history policy adapted from deep_autolabel memory_mode/spatial_sliding.
#include "LIVMapper.h"
#include "memory_mode/page/voxel_pager.h"
#include <algorithm>

size_t VoxelMapManager::ApplySpatialWindow(const flivo::SpatialWindowBounds& bounds) {
  spatial_window_=bounds;
  size_t removed=0;
  for(auto it=voxel_map_.begin();it!=voxel_map_.end();) {
    const auto& key=it->first;
    if(bounds.IntersectsVoxel(key.x,key.y,key.z,config_setting_.max_voxel_size_)) {++it;continue;}
    if(pager_)pager_->Forget(key);
    delete it->second;it=voxel_map_.erase(it);++removed;
  }
  ptpl_list_.clear();
  return removed;
}

size_t VIOManager::ApplySpatialWindow(const flivo::SpatialWindowBounds& bounds) {
  for(auto& camera:cameras)camera->spatial_window=bounds;
  const auto outside=[&](const auto& item) {
    const auto& key=item.first;
    return !bounds.IntersectsVoxel(key.x,key.y,key.z,.5);
  };
  if(std::none_of(visual_map.begin(),visual_map.end(),outside))return 0;
  // All cameras borrow points from the same map. Clear every camera before
  // destroying a single owning voxel, after all solver/export work has ended.
  for(auto& camera:cameras)camera->clearStartupMap();
  size_t removed=0;
  for(auto it=visual_map.begin();it!=visual_map.end();) {
    if(!outside(*it)) {++it;continue;}
    removed+=it->second->voxel_points.size();
    delete it->second;it=visual_map.erase(it);
  }
  return removed;
}

void LIVMapper::finishMemoryFrame() {
  if(!spatial_window || !lidar_map_inited || (imu_en && p_imu->imu_need_init) ||
      startup_replay.pending() || startup_replay.prepared() || startup_replay.replaying())return;
  // Image-clock LIVO has an LIO half followed by its VIO/export half. A scan
  // clock or LiDAR-only input instead completes each independently timed event.
  if(img_en && lidar_en && !scan_clock && LidarMeasures.lio_vio_flg!=VIO)return;
  if(!spatial_window->Update(_state.pos_end))return;
  const auto& bounds=*spatial_window->bounds();
  sliding_vio_evictions+=vio_manager->ApplySpatialWindow(bounds);
  sliding_lio_evictions+=voxelmap_manager->ApplySpatialWindow(bounds);
  ++window_moves;
}
