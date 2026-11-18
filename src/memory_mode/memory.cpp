#include "LIVMapper.h"
#include <sys/resource.h>

void LIVMapper::configureMemory() {
  if(memory_mode=="spatial_sliding")
    spatial_window=std::make_unique<flivo::SpatialSlidingWindow>(sliding_half_extent_m,sliding_trigger_m);
  if(memory_mode=="page") {
    voxelmap_manager->EnablePaging(DEBUG_FILE_DIR(".memory_cache/lio.bin"));
    for(auto& camera:vio_manager->cameras)
      camera->plane_lookup=[this](const VOXEL_LOCATION& key){return voxelmap_manager->GetVoxel(key);};
  }
  // Spatial pruning bounds geometry, while retained reference images can still
  // fill RAM with several cameras. Back the same pixels without changing which
  // observations or voxels survive the shared spatial window.
  if(memory_mode=="page" || memory_mode=="spatial_sliding") {
    for(auto& camera:vio_manager->cameras) {
      camera->image_store=std::make_unique<flivo::ImageBackingStore>(
          std::filesystem::path(DEBUG_FILE_DIR(".memory_cache/images"))/camera->id);
    }
  }
  memory_log.open(DEBUG_FILE_DIR("memory_working_set.csv"));
  memory_log<<"time_s,mode,rss_kib,peak_rss_kib,lio_total,lio_resident,lio_query_working,lio_loads,lio_evictions,lio_bytes_written,vio_voxels,vio_points,vio_observations,vio_reference_images,vio_mapped_image_bytes,vio_working_image_bytes,vio_image_drop_calls,window_moves,sliding_lio_evictions,sliding_vio_evictions,sliding_lio_rejected_points,sliding_vio_rejected_points,window_active,window_x,window_y,window_z,window_half_extent_m\n";
  if(!memory_log)throw std::runtime_error("cannot create memory working-set log");
}

void LIVMapper::recordMemory(bool force) {
  const double time=LidarMeasures.last_lio_update_time;
  if(!memory_log.is_open() || (!force && memory_last_recorded>=0 && time-memory_last_recorded<1.))return;
  memory_last_recorded=time;
  size_t points=0,observations=0,mapped=0,working=0;uint64_t drops=0,rejected=0;
  std::set<const unsigned char*> images;
  for(const auto& voxel:vio_manager->visual_map)for(const auto* point:voxel.second->voxel_points) {
    ++points;
    for(const auto& camera:point->camera_observations)for(const auto* observation:camera.features) {
      ++observations;images.insert(observation->img_.data);
    }
  }
  for(const auto& camera:vio_manager->cameras)if(camera->image_store) {
    const auto& s=camera->image_store->stats;
    mapped+=s.mapped_bytes;working+=s.working_bytes;drops+=s.drop_calls;
  }
  for(const auto& camera:vio_manager->cameras)rejected+=camera->sliding_rejected_points;
  const bool active=spatial_window && spatial_window->bounds().has_value();
  const V3D center=active ? spatial_window->bounds()->center : V3D::Zero();
  size_t pages=0,resident=0;std::ifstream stat("/proc/self/statm");
  if(!(stat>>pages>>resident))throw std::runtime_error("cannot read process RSS");
  struct rusage usage{};
  if(getrusage(RUSAGE_SELF,&usage)!=0)throw std::runtime_error("cannot read process peak RSS");
  const auto c=voxelmap_manager->GetCacheStats();
  memory_log<<std::setprecision(17)<<time<<','<<memory_mode<<','
    <<resident*static_cast<size_t>(sysconf(_SC_PAGESIZE))/1024<<','<<usage.ru_maxrss
    <<','<<c.total<<','<<c.resident<<','<<c.working<<','<<c.loads<<','<<c.evictions<<','<<c.bytes_written
    <<','<<vio_manager->visual_map.size()<<','<<points<<','<<observations<<','<<images.size()
    <<','<<mapped<<','<<working<<','<<drops
    <<','<<window_moves<<','<<sliding_lio_evictions<<','<<sliding_vio_evictions
    <<','<<voxelmap_manager->sliding_rejected_points_<<','<<rejected<<','<<active
    <<','<<center.x()<<','<<center.y()<<','<<center.z()<<','<<sliding_half_extent_m<<'\n';
  memory_log.flush();
  if(!memory_log)throw std::runtime_error("cannot write memory working-set log");
}
