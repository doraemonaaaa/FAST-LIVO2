#include "LIVMapper.h"
#include "imu_support.h"
#include <iomanip>

void LIVMapper::consumeLidarFront() {
  last_lidar_support_begin=lid_header_time_buffer.front();
  lid_raw_data_buffer.pop_front();lid_header_time_buffer.pop_front();
  if(imu_checked_scans)--imu_checked_scans;
}

void LIVMapper::filterImuCoveredScans() {
  bool eof=false;ros::param::get("/fast_livo/input_finished",eof);
  while(imu_checked_scans<lid_raw_data_buffer.size()) {
    const size_t i=imu_checked_scans;auto& cloud=*lid_raw_data_buffer[i];
    if(!std::is_sorted(cloud.points.begin(),cloud.points.end(),[](const auto& a,const auto& b){return a.curvature<b.curvature;}))
      std::stable_sort(cloud.points.begin(),cloud.points.end(),[](const auto& a,const auto& b){return a.curvature<b.curvature;});
    const double begin=lid_header_time_buffer[i];double end=begin;
    for(const auto& p:cloud)end=std::max(end,begin+p.curvature*.001);
    const auto status=cloud.empty()?imu_support::Coverage::Invalid:
      imu_support::Check(scan_imu_history,begin,end);
    if(status==imu_support::Coverage::MissingEnd && !eof)break;
    if(imu_coverage_log) {
      imu_coverage_log<<std::setprecision(17)<<begin<<','<<end<<','<<cloud.size()<<','
        <<(status==imu_support::Coverage::Covered?"accepted":"dropped")<<','<<imu_support::Reason(status)<<',';
      if(!scan_imu_history.empty())imu_coverage_log<<imu_support::Time(scan_imu_history.front())<<','<<imu_support::Time(scan_imu_history.back());
      else imu_coverage_log<<',';
      imu_coverage_log<<'\n'<<std::flush;
    }
    if(status==imu_support::Coverage::Covered) {
      startup_replay.observeScan(cloud,begin);++imu_checked_scans;continue;
    }
    ROS_WARN("[IMU coverage] discard entire scan %.9f..%.9f (%zu points): %s",begin,end,cloud.size(),imu_support::Reason(status));
    lid_raw_data_buffer.erase(lid_raw_data_buffer.begin()+i);
    lid_header_time_buffer.erase(lid_header_time_buffer.begin()+i);
  }
  if(imu_checked_scans && LidarMeasures.last_lio_update_time<0)
    LidarMeasures.last_lio_update_time=lid_header_time_buffer.front();
  // Keep a real left bracket for both the next scan and the shared estimator.
  double retain=p_imu->propagationTime();
  if(last_lidar_support_begin>=0)retain=std::min(retain,last_lidar_support_begin);
  if(!lid_header_time_buffer.empty())retain=std::min(retain,lid_header_time_buffer.front());
  while(scan_imu_history.size()>2 && imu_support::Time(scan_imu_history[1])<=retain)
    scan_imu_history.pop_front();
}

deque<sensor_msgs::Imu::ConstPtr> LIVMapper::imuBracket(double begin,double end) const {
  if(imu_support::Check(scan_imu_history,begin,end,std::numeric_limits<double>::infinity())!=imu_support::Coverage::Covered)
    throw std::runtime_error("State prediction lacks real IMU endpoint brackets");
  auto first=std::upper_bound(scan_imu_history.begin(),scan_imu_history.end(),begin,
    [](double t,const auto& s){return t<imu_support::Time(s);});
  --first;
  deque<sensor_msgs::Imu::ConstPtr> out;
  for(auto p=first;p!=scan_imu_history.end();++p) {
    out.push_back(*p);if(imu_support::Time(*p)>=end)break;
  }
  return out;
}

bool LIVMapper::sync_packages(LidarMeasureGroup& meas) {
  if(imu_en)filterImuCoveredScans();
  if(!syncRawPackages(meas))return false;
  if(imu_en) {
    auto& m=meas.measures.back();const double end=meas.lio_vio_flg==LIO?m.lio_time:m.vio_time;
    const double begin=last_sync_time<0?meas.last_lio_update_time:
      std::min(meas.last_lio_update_time,last_sync_time);
    // Peek at the first real IMU at/after the target; it remains available for
    // the next event. Within-range interpolation replaces the former hold.
    m.imu=imuBracket(begin,end);
    last_sync_time=end;
  }
  return true;
}
