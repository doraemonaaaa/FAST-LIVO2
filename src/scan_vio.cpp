#include "LIVMapper.h"
#include "camera_events.h"
#include <limits>
#include <iomanip>
#include <stdexcept>

bool LIVMapper::syncScanEvent(LidarMeasureGroup &meas)
{
  bool eof=false;ros::param::get("/fast_livo/input_finished",eof);
  if(meas.last_lio_update_time<0) {
    if(lid_header_time_buffer.empty())return false;
    meas.last_lio_update_time=lid_header_time_buffer.front();
  }
  double end=std::numeric_limits<double>::infinity();
  if(!lid_raw_data_buffer.empty()) {
    scan_begin=lid_header_time_buffer.front();
    end=scan_begin+lid_raw_data_buffer.front()->points.back().curvature*.001;
    camera_lidar_coverage=std::max(camera_lidar_coverage,
      lid_header_time_buffer.back()+lid_raw_data_buffer.back()->points.back().curvature*.001);
  }
  const double coverage=std::min(camera_lidar_coverage,last_timestamp_imu);
  double image_time;
  const bool image_ready=prepareCameraEvent(image_time,coverage);
  bool scan_ready=end<=coverage;
  if(scan_ready) {
    std::lock_guard<std::mutex> lock(mtx_buffer);
    vector<double> first,last;
    for(const auto &q:camera_queues) {
      first.push_back(q.times.empty()?std::numeric_limits<double>::infinity():q.times.front());
      last.push_back(q.last_timestamp);
    }
    scan_ready=cameraEventReady(end,first,last,coverage,camera_maximum_lateness,eof);
    if(std::isfinite(image_time) && image_time<end && !image_ready)scan_ready=false;
  }
  if(!image_ready && !scan_ready)return false;
  MeasureGroup m;
  bundle_sequential_updates=false;
  meas.pcl_proc_cur->clear();
  if(scan_ready && (!image_ready || end<=image_time)) {
    std::lock_guard<std::mutex> lock(mtx_buffer);
    *meas.pcl_proc_cur=*lid_raw_data_buffer.front();
    lid_raw_data_buffer.pop_front();lid_header_time_buffer.pop_front();
    meas.last_lio_update_time=scan_begin; // Raw point offsets use this scan's origin.
    m.lio_time=end;m.vio_time=end;meas.lio_vio_flg=LIO;
    sequential_anchor_event=false;pending_images.clear();
  } else {
    m.lio_time=meas.last_lio_update_time;m.vio_time=image_time;
    m.images=pending_images;m.img=pending_images.front().pixels;
    meas.lio_vio_flg=VIO;sequential_anchor_event=false;
    {std::lock_guard<std::mutex> lock(mtx_buffer);consumeImageBundle();}
    last_camera_event_stamp=m.images.front().timestamp_ns;
    last_camera_event_time=image_time;
  }
  const double target=meas.lio_vio_flg==LIO?m.lio_time:m.vio_time;
  {
    std::lock_guard<std::mutex> lock(mtx_buffer);
    while(!imu_buffer.empty() && imu_buffer.front()->header.stamp.toSec()<=target) {
      m.imu.push_back(imu_buffer.front());imu_buffer.pop_front();
    }
  }
  meas.measures.clear();meas.measures.push_back(m);sig_buffer.notify_all();return true;
}

void LIVMapper::processScanImu()
{
  const bool lidar=LidarMeasures.lio_vio_flg==LIO;
  if(p_imu->imu_need_init) {
    p_imu->Process2(LidarMeasures,_state,feats_undistort);
  } else {
    // Propagate the one shared estimator. No raw scan points are partitioned at images.
    LidarMeasureGroup propagation=LidarMeasures;
    propagation.pcl_proc_cur.reset(new PointCloudXYZI());
    p_imu->Process2(propagation,_state,feats_undistort);
    LidarMeasures.last_lio_update_time=propagation.last_lio_update_time;
    if(lidar) {
      deque<sensor_msgs::Imu::ConstPtr> samples;
      const double end=LidarMeasures.measures.back().lio_time;
      {
        std::lock_guard<std::mutex> lock(mtx_buffer);
        while(scan_imu_history.size()>1 && scan_imu_history[1]->header.stamp.toSec()<=scan_begin)
          scan_imu_history.pop_front();
        for(const auto &sample:scan_imu_history)
          if(sample->header.stamp.toSec()<=end)samples.push_back(sample);
      }
      scan_audit<<std::setprecision(17)<<scan_begin<<','<<end<<','<<LidarMeasures.pcl_proc_cur->size()<<','<<samples.size()<<','
        <<_state.vel_end.norm()<<','<<_state.bias_g.norm()<<','<<_state.bias_a.norm()<<','<<_state.gravity.norm()<<'\n'<<std::flush;
      p_imu->deskewWholeScan(*LidarMeasures.pcl_proc_cur,scan_begin,end,_state,samples,feats_undistort);
      if(feats_undistort->size()!=LidarMeasures.pcl_proc_cur->size())
        throw std::runtime_error("Whole-scan deskew changed point count");
    }
  }
  if(gravity_align_en)gravityAlignment();
  state_propagat=_state;voxelmap_manager->state_=_state;
  voxelmap_manager->feats_undistort_=feats_undistort;
}
