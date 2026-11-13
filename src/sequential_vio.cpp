#include "LIVMapper.h"
#include "imu_support.h"
#include "livo_time_cut.h"
#include "camera_events.h"
#include "pending_cloud.h"
#include <limits>
#include <iomanip>
#include <stdexcept>

bool LIVMapper::prepareCameraEvent(double &timestamp,double sensor_time)
{
  std::lock_guard<std::mutex> lock(mtx_buffer);
  bool eof=false;ros::param::get("/fast_livo/input_finished",eof);
  vector<double> first,last;
  size_t selected=camera_queues.size();timestamp=std::numeric_limits<double>::infinity();
  for(size_t k=0;k<camera_queues.size();++k) {
    auto &q=camera_queues[k];
    while(!q.times.empty() && (q.times.front()<LidarMeasures.last_lio_update_time || q.stamps.front()<last_camera_event_stamp ||
          (eof && q.times.front()>sensor_time) ||
          (q.times.front()<=last_timestamp_imu && imu_support::Check(scan_imu_history,q.times.front(),q.times.front())!=imu_support::Coverage::Covered))) {
      camera_input_log<<vio_manager->cameras[k]->id<<','<<q.stamps.front()<<','
          <<(q.times.front()<LidarMeasures.last_lio_update_time || q.stamps.front()<last_camera_event_stamp?"before_state":"outside_sensor_coverage")<<'\n'<<std::flush;
      q.times.pop_front();q.images.pop_front();q.stamps.pop_front();
    }
    const double t=q.times.empty()?std::numeric_limits<double>::infinity():q.times.front();
    first.push_back(t);last.push_back(q.last_timestamp);
    if(std::isfinite(t) && (selected==camera_queues.size() || q.stamps.front()<camera_queues[selected].stamps.front()))
      {timestamp=t;selected=k;}
  }
  if(selected==camera_queues.size() || timestamp>sensor_time ||
      !cameraEventReady(timestamp,first,last,sensor_time,camera_maximum_lateness,eof))return false;
  auto &q=camera_queues[selected];
  pending_images={{selected,q.times.front(),q.images.front(),q.stamps.front()}};
  sequential_anchor_event=selected==vio_manager->anchor_index;
  return true;
}

bool LIVMapper::syncCameraEvent(LidarMeasureGroup &meas)
{
  if(meas.lio_vio_flg==LIO) {
    // The anchor event was reserved before its LIO. No fresh input is needed.
    MeasureGroup m;m.lio_time=meas.last_lio_update_time;m.vio_time=pending_images.front().timestamp;
    m.img=pending_images.front().pixels;m.images=pending_images;
    meas.measures.clear();meas.measures.push_back(m);meas.lio_vio_flg=VIO;
    {std::lock_guard<std::mutex> lock(mtx_buffer);consumeImageBundle();}
    last_camera_event_time=m.vio_time;last_camera_event_stamp=m.images.front().timestamp_ns;
    sig_buffer.notify_all();return true;
  }
  if(meas.last_lio_update_time<0) {
    if(lid_header_time_buffer.empty() || !imu_checked_scans)return false;
    meas.last_lio_update_time=lid_header_time_buffer.front();
  }
  if(imu_checked_scans)camera_lidar_coverage=std::max(camera_lidar_coverage,
      lid_header_time_buffer[imu_checked_scans-1]+lid_raw_data_buffer[imu_checked_scans-1]->points.back().curvature*.001);
  const double coverage=std::min(camera_lidar_coverage,last_timestamp_imu);
  double timestamp;
  if(!prepareCameraEvent(timestamp,coverage))return false;
  const double begin=meas.last_lio_update_time;
  MeasureGroup m=collectLivoSlice(meas,timestamp);
  meas.measures.clear();
  bundle_sequential_updates=!sequential_anchor_event;
  if(sequential_anchor_event) {
    bundle_lio_time=timestamp;meas.lio_vio_flg=LIO;
  } else {
    m.lio_time=begin;m.vio_time=timestamp;m.img=pending_images.front().pixels;m.images=pending_images;
    meas.lio_vio_flg=VIO;
    {std::lock_guard<std::mutex> lock(mtx_buffer);consumeImageBundle();}
    last_camera_event_time=timestamp;
    last_camera_event_stamp=m.images.front().timestamp_ns;
  }
  meas.measures.push_back(m);sig_buffer.notify_all();return true;
}

MeasureGroup LIVMapper::collectLivoSlice(LidarMeasureGroup &meas,double img_capture_time)
{
  struct MeasureGroup m;

  // printf("[ Data Cut ] LIO \n");
  // printf("[ Data Cut ] img_capture_time: %lf \n", img_capture_time);
  m.imu.clear();
  m.lio_time = img_capture_time;
  mtx_buffer.lock();
  while (!imu_buffer.empty())
  {
    if (imu_buffer.front()->header.stamp.toSec() > m.lio_time) break;

    if (imu_buffer.front()->header.stamp.toSec() > meas.last_lio_update_time) m.imu.push_back(imu_buffer.front());

    imu_buffer.pop_front();
    // printf("[ Data Cut ] imu time: %lf \n",
    // imu_buffer.front()->header.stamp.toSec());
  }
  mtx_buffer.unlock();
  sig_buffer.notify_all();

  splitPendingPointsAtImageTime(meas.pcl_proc_next->points, meas.pcl_proc_cur->points,
                               (img_capture_time - meas.last_lio_update_time) * 1000.0);

  int lid_frame_num = lid_raw_data_buffer.size();
  int max_size = meas.pcl_proc_cur->size() + 24000 * lid_frame_num;
  meas.pcl_proc_cur->reserve(max_size);
  meas.pcl_proc_next->reserve(max_size);
  // deque<PointCloudXYZI::Ptr> lidar_buffer_tmp;

  while (!lid_raw_data_buffer.empty() && (!imu_en || imu_checked_scans))
  {
    if (lid_header_time_buffer.front() > img_capture_time) break;
    auto pcl(lid_raw_data_buffer.front()->points);
    double frame_header_time(lid_header_time_buffer.front());
    float max_offs_time_ms = (m.lio_time - frame_header_time) * 1000.0f;

    for (int i = 0; i < pcl.size(); i++)
    {
      auto pt = pcl[i];
      if (pcl[i].curvature < max_offs_time_ms)
      {
        pt.curvature += (frame_header_time - meas.last_lio_update_time) * 1000.0f;
        meas.pcl_proc_cur->points.push_back(pt);
      }
      else
      {
        pt.curvature += (frame_header_time - m.lio_time) * 1000.0f;
        meas.pcl_proc_next->points.push_back(pt);
      }
    }
    consumeLidarFront();
  }

  // Individually sorted scans can overlap in time. Deskew requires the
  // combined current cloud to be sorted as well.
  std::stable_sort(meas.pcl_proc_cur->points.begin(), meas.pcl_proc_cur->points.end(),
                   [](const PointType &a, const PointType &b) { return a.curvature < b.curvature; });
  std::stable_sort(meas.pcl_proc_next->points.begin(), meas.pcl_proc_next->points.end(),
                   [](const PointType &a, const PointType &b) { return a.curvature < b.curvature; });

  return m;
}

void LIVMapper::processSequentialVIO(bool visual_support)
{
  const auto &measurement=LidarMeasures.measures.back();
  const double origin=measurement.lio_time;
  double cursor=origin;
  size_t imu_index=0,point_index=0;
  const auto &points=LidarMeasures.pcl_proc_cur->points;
  M3D visual_prior_rotation;
  V3D visual_prior_position;
  auto advance=[&](double timestamp) {
    if(timestamp<cursor)throw std::runtime_error("Non-monotonic visual state propagation");
    LidarMeasureGroup slice;slice.last_lio_update_time=cursor;
    // During the preparation pass only LIO slices may run the motion fit.
    slice.lio_vio_flg=p_imu->imu_need_init?VIO:LIO;
    MeasureGroup m;m.lio_time=timestamp;m.vio_time=timestamp;
    while(imu_index<measurement.imu.size() && measurement.imu[imu_index]->header.stamp.toSec()<=timestamp)
      m.imu.push_back(measurement.imu[imu_index++]);
    if(imu_index<measurement.imu.size())m.imu.push_back(measurement.imu[imu_index]);
    vector<double> point_times;
    while(point_index<points.size() &&
          (timestamp==measurement.vio_time || points[point_index].curvature*.001<=timestamp-origin)) {
      PointType point=points[point_index++];
      const double offset=point.curvature*.001,time=origin+offset;
      point.curvature=(offset-(cursor-origin))*1000.;
      slice.pcl_proc_cur->push_back(point);point_times.push_back(time);
    }
    slice.measures.push_back(m);
    PointCloudXYZI::Ptr deskewed(new PointCloudXYZI());
    constrainGravityBeforePrediction(timestamp);
    p_imu->Process2(slice,_state,deskewed,true);
    if(!p_imu->imu_need_init && !startup_replay.pending()) {
      if(deskewed->size()!=point_times.size())throw std::runtime_error("Inter-camera LiDAR point count changed");
      for(size_t i=0;i<deskewed->size();++i) {
        PointType point=deskewed->points[i];
        const V3D world=_state.rot_end*(extR*V3D(point.x,point.y,point.z)+extT)+_state.pos_end;
        inter_camera_points.push_back({point,world,point_times[i]});
      }
    }
    state_propagat=_state;
    visual_prior_rotation=_state.rot_end;
    visual_prior_position=_state.pos_end;
    LidarMeasures.last_lio_update_time=timestamp;
    camera_events_log<<std::setprecision(17)<<cursor<<','<<timestamp<<','<<m.imu.size()<<','
      <<slice.pcl_proc_cur->size()<<','<<!p_imu->imu_need_init<<','<<_state.pos_end.x()<<','
      <<_state.pos_end.y()<<','<<_state.pos_end.z()<<','<<_state.cov.trace()<<'\n'<<std::flush;
    cursor=timestamp;
  };
  if(visual_support) {
    // The independent scheduler reserves exactly one acquisition per call.
    // Correct the entire not-yet-matched scan, including earlier camera slices.
    if(measurement.images.size()!=1)throw std::runtime_error("Expected one independent camera event");
    vio_manager->processFrame(measurement.images,_pv_list,voxelmap_manager->voxel_map_,origin,advance);
    correctPendingCloud(inter_camera_points,visual_prior_rotation,visual_prior_position,
                        _state.rot_end,_state.pos_end);
  }
  else
    for(size_t i:cameraEventOrder(measurement.images))advance(measurement.images[i].timestamp);
  if(point_index!=points.size() || (imu_index<measurement.imu.size() && measurement.imu[imu_index]->header.stamp.toSec()<=measurement.vio_time))
    throw std::runtime_error("Sequential bundle left unprocessed sensor input");
}
