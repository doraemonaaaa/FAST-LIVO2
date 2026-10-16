#include "replay/startup_replay.h"
#include "replay/initialization_checkpoint.h"
#include "LIVMapper.h"
#include <map>
#include <cstdio>
#include <ros/serialization.h>

void StartupReplay::readParameters(ros::NodeHandle &nh, bool imu_en)
{
  bool replay_enabled;
  nh.param<bool>("imu/initialization_replay",replay_enabled,false);
  if (!replay_enabled) return;
  int maximum_mb;
  std::string phase;
  nh.param<int>("imu/initialization_replay_max_mb",maximum_mb,2048);
  nh.param<std::string>("imu/initialization_replay_phase",phase,"");
  nh.param<std::string>("imu/initialization_checkpoint",checkpoint_,"");
  nh.param<std::string>("imu/initialization_identity",identity_,"");
  if (!imu_en || maximum_mb<=0 || (phase!="prepare" && phase!="run") || checkpoint_.empty() || identity_.empty())
    throw std::runtime_error("Two-pass initialization requires the offline export runner, enabled IMU and a positive cache budget");
  phase_=phase=="prepare" ? Phase::Preparing : Phase::Replaying;
  limit_=size_t(maximum_mb)*1024*1024;
}

void StartupReplay::initializeLog(const std::string &path)
{
  if (!enabled()) return;
  log_.open(path);
  if (!log_) throw std::runtime_error("Cannot create initialization replay log");
  log_ << "event,begin,end,groups,bytes\n" << std::flush;
  input_log_.open(path.substr(0,path.find_last_of('/'))+"/initialization_inputs.csv");
  if (!input_log_) throw std::runtime_error("Cannot create initialization input audit");
  input_log_ << "index,begin,end,flag,points,images,fingerprint\n" << std::flush;
}

void StartupReplay::retainImuForInitialization(const sensor_msgs::Imu::ConstPtr &imu)
{
  if (!pending()) return;
  size_t bytes=sizeof(sensor_msgs::Imu)+imu->header.frame_id.size();
  if (bytes>limit_ || bytes_>limit_-bytes) throw std::runtime_error("Initialization IMU/index budget exceeded");
  imu_samples_.push_back(imu); bytes_+=bytes;
}

StartupReplay::Index StartupReplay::index(const LidarMeasureGroup &group) const
{
  const auto &m=group.measures.back();
  Index result{group.last_lio_update_time,group.lio_vio_flg==LIO?m.lio_time:m.vio_time,
      1469598103934665603ULL,group.pcl_proc_cur->size(),m.images.size(),int(group.lio_vio_flg)};
  auto append=[&](const void *p,size_t n){result.fingerprint=flivo_trace::hashBytes(p,n,result.fingerprint);};
  append(&result.flag,sizeof(result.flag)); append(&result.begin,sizeof(result.begin)); append(&result.end,sizeof(result.end));
  for (const auto &p : group.pcl_proc_cur->points)
  { const float values[]={p.x,p.y,p.z,p.intensity,p.curvature}; append(values,sizeof(values)); }
  for (const auto &imu : m.imu)
  {
    std::vector<uint8_t> data(ros::serialization::serializationLength(*imu));
    ros::serialization::OStream stream(data.data(),data.size()); ros::serialization::serialize(stream,*imu);
    append(data.data(),data.size());
  }
  for (const auto &image : m.images)
  {
    append(&image.camera_index,sizeof(image.camera_index)); append(&image.timestamp_ns,sizeof(image.timestamp_ns));
    append(&image.timestamp,sizeof(image.timestamp));
    const int shape[]={image.pixels.rows,image.pixels.cols,image.pixels.type()}; append(shape,sizeof(shape));
    for (int row=0;row<image.pixels.rows;++row) append(image.pixels.ptr(row),image.pixels.cols*image.pixels.elemSize());
  }
  return result;
}

void StartupReplay::observeMeasurement(const LidarMeasureGroup &group)
{
  if (!pending() && !replaying()) return;
  const Index current=index(group);
  if (pending())
  {
    if (indices_.empty()) begin_=current.begin;
    if (bytes_>limit_-sizeof(Index)) throw std::runtime_error("Initialization index budget exceeded");
    indices_.push_back(current); bytes_+=sizeof(Index);
    log_ << std::setprecision(17) << "index," << current.begin << ',' << current.end << ','
         << indices_.size() << ',' << bytes_ << '\n' << std::flush;
  }
  else
  {
    const auto &expected=indices_.at(verified_);
    if (current.begin!=expected.begin || current.end!=expected.end || current.flag!=expected.flag ||
        current.points!=expected.points || current.images!=expected.images || current.fingerprint!=expected.fingerprint)
    {
      std::ostringstream detail;
      detail << std::setprecision(17) << "Initialization reread input mismatch at group " << verified_+1
             << "; origin expected=" << expected.begin << " actual=" << current.begin
             << "; end expected=" << expected.end << " actual=" << current.end
             << "; points expected=" << expected.points << " actual=" << current.points
             << "; images expected=" << expected.images << " actual=" << current.images
             << "; hash expected=" << expected.fingerprint << " actual=" << current.fingerprint;
      throw std::runtime_error(detail.str());
    }
  }
  input_log_ << std::setprecision(17) << (pending()?indices_.size():verified_+1) << ','
             << current.begin << ',' << current.end << ',' << current.flag << ','
             << current.points << ',' << current.images << ',' << current.fingerprint << '\n' << std::flush;
}

void StartupReplay::prepareNextInput(LidarMeasureGroup &group) const
{
  // Initialization can leave the synchronization cursor unchanged when the
  // leading slice has no IMU. The restored estimator already integrates such
  // a slice, so reproduce the recorded slicing origin independently of its
  // updated state. This changes no IMU integration history or EKF state.
  if (replaying())
    group.last_lio_update_time=indices_[verified_].begin;
}

void StartupReplay::afterUpdate(LIVMapper &mapper)
{
  if (replaying())
  {
    ++verified_;
    log_ << std::setprecision(17) << "group," << mapper.LidarMeasures.lidar_frame_beg_time << ','
         << mapper.LidarMeasures.last_lio_update_time << ',' << verified_ << ",0\n";
    if (verified_==indices_.size())
    {
      log_ << "complete," << begin_ << ',' << end_ << ',' << indices_.size() << ',' << bytes_ << '\n' << std::flush;
      phase_=Phase::Complete;
    }
    return;
  }
  if (!pending() || mapper.p_imu->imu_need_init) return;
  end_=mapper.LidarMeasures.last_lio_update_time;
  std::map<uint64_t,sensor_msgs::Imu::ConstPtr> ordered;
  for (const auto &imu : imu_samples_)
    if (imu->header.stamp.toSec()<=end_) ordered.emplace(imu->header.stamp.toNSec(),imu);
  std::deque<sensor_msgs::Imu::ConstPtr> samples;
  for (const auto &item : ordered)
  { if (item.second->header.stamp.toSec()<=begin_) samples.clear(); samples.push_back(item.second); }
  mapper.p_imu->rewindForStartupReplay(mapper._state,end_,begin_,samples);
  std::ostringstream out;
  out << std::setprecision(17) << std::quoted(identity_) << '\n' << begin_ << ' ' << end_ << ' '
      << indices_.size() << ' ' << bytes_ << ' ' << mapper.gravity_align_finished << '\n';
  initialization_checkpoint::writeState(out,mapper._state);
  mapper.p_imu->writeInitialization(out);
  for (const auto &i : indices_) out << i.begin << ' ' << i.end << ' ' << i.fingerprint << ' ' << i.points << ' ' << i.images << ' ' << i.flag << '\n';
  const auto payload=out.str(); const auto temporary=checkpoint_+".tmp";
  { std::ofstream file(temporary,std::ios::binary);
    file << "FLIVO_INITIALIZATION_V1\n" << payload.size() << ' ' << initialization_checkpoint::hash(payload) << '\n' << payload;
    file.flush(); if (!file) throw std::runtime_error("Cannot write initialization checkpoint"); }
  if (std::rename(temporary.c_str(),checkpoint_.c_str())) throw std::runtime_error("Cannot commit initialization checkpoint");
  log_ << "prepared," << begin_ << ',' << end_ << ',' << indices_.size() << ',' << bytes_ << '\n' << std::flush;
  phase_=Phase::Prepared;
  ROS_WARN("[initialization] checkpoint ready: %zu indices, %zu IMU samples, %.3f MiB retained, no point/image payload cache",
      indices_.size(),imu_samples_.size(),bytes_/1048576.);
  imu_samples_.clear();
}

void StartupReplay::restore(LIVMapper &mapper)
{
  if (!replaying()) return;
  std::ifstream file(checkpoint_,std::ios::binary); std::string magic; size_t length=0; uint64_t checksum=0;
  std::getline(file,magic);
  if (magic!="FLIVO_INITIALIZATION_V1" || !(file >> length >> checksum) || length==0 || length>64*1024*1024 || file.get()!='\n')
    throw std::runtime_error("Invalid initialization checkpoint header");
  std::string payload(length,'\0');file.read(&payload[0],length);
  if (size_t(file.gcount())!=length || file.peek()!=EOF || initialization_checkpoint::hash(payload)!=checksum)
    throw std::runtime_error("Truncated or corrupt initialization checkpoint");
  std::istringstream in(payload);std::string identity; size_t count=0;
  if (!(in >> std::quoted(identity)) || identity!=identity_ ||
      !(in >> begin_ >> end_ >> count >> bytes_ >> mapper.gravity_align_finished) ||
      !std::isfinite(begin_) || !std::isfinite(end_) || end_<begin_ || count==0 || count>1000000)
    throw std::runtime_error("Initialization checkpoint input/config/runtime identity mismatch");
  initialization_checkpoint::readState(in,mapper._state);
  mapper.p_imu->readInitialization(in);
  indices_.resize(count);
  for (auto &i : indices_)
    if (!(in >> i.begin >> i.end >> i.fingerprint >> i.points >> i.images >> i.flag))
      throw std::runtime_error("Invalid initialization input index");
  in >> std::ws;if (!in.eof()) throw std::runtime_error("Trailing initialization checkpoint data");
  mapper.state_propagat=mapper._state; mapper.voxelmap_manager->state_=mapper._state;
  log_ << std::setprecision(17) << "begin," << begin_ << ',' << end_ << ',' << indices_.size() << ',' << bytes_ << '\n' << std::flush;
  ROS_WARN("[initialization] restored start state at %.9f; verifying %zu reread input groups",begin_,indices_.size());
}

void StartupReplay::finish()
{
  log_.flush();
  input_log_.flush();
  if (pending() || replaying())
    throw std::runtime_error("Initialization preparation/reread did not complete");
}
