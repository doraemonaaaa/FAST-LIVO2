#include "LIVMapper.h"
#include "camera_sync.h"
#include <vikit/camera_loader.h>
#include <regex>
#include <filesystem>
#include <iomanip>
#include <set>
#include <stdexcept>

void LIVMapper::configureCameras()
{
  ros::NodeHandle nh;
  XmlRpc::XmlRpcValue definitions;
  vector<string> order;
  const bool supplied = nh.getParam("cameras", definitions);
  if (supplied)
  {
    if (definitions.getType() != XmlRpc::XmlRpcValue::TypeStruct || definitions.size() == 0)
      throw std::runtime_error("cameras must be a nonempty map");
    if (!nh.getParam("camera_order", order))
      for (auto it = definitions.begin(); it != definitions.end(); ++it) order.push_back(it->first);
  }
  else order.push_back("camera0"); // Normalize legacy configuration, not an algorithm branch.
  if (order.empty()) throw std::runtime_error("camera_order must not be empty");
  string anchor;
  nh.param<string>("anchor_camera", anchor, order.front());
  nh.param<double>("camera_sync/tolerance_s", vio_manager->sync_tolerance, 0.0001);
  nh.param<double>("camera_sync/maximum_lateness_s", camera_maximum_lateness, 0.5);
  if (!(vio_manager->sync_tolerance >= 0 && vio_manager->sync_tolerance <= .005 &&
        camera_maximum_lateness > 2 * vio_manager->sync_tolerance && camera_maximum_lateness < 2))
    throw std::runtime_error("Invalid synchronous camera tolerance/lateness");
  std::set<string> ids, topics;
  for (size_t k = 0; k < order.size(); ++k)
  {
    const string &id = order[k];
    if (!std::regex_match(id, std::regex("[A-Za-z0-9_]+")) || !ids.insert(id).second)
      throw std::runtime_error("Invalid/duplicate camera id");
    if (supplied && !definitions.hasMember(id)) throw std::runtime_error("Camera not defined: " + id);
    const string prefix = supplied ? "cameras/" + id + "/" : "";
    auto view = std::make_unique<VisualCamera>(vio_manager->visual_map);
    view->id = id; view->camera_index = k; view->camera_count = order.size();
    if (supplied)
    {
      if (!nh.getParam(prefix + "topic", view->topic)) throw std::runtime_error("Missing camera topic: " + id);
    }
    else view->topic = img_topic;
    if (view->topic.empty() || view->topic[0] != '/' || !topics.insert(view->topic).second)
      throw std::runtime_error("Invalid/duplicate camera topic");
    string camera_ns = supplied ? "laserMapping/cameras/" + id : "laserMapping";
    if (!vk::camera_loader::loadFromRosNs(camera_ns, view->cam))
      throw std::runtime_error("Camera model not correctly specified: " + id);
    vector<double> rotation = cameraextrinR, translation = cameraextrinT;
    if (supplied && (!nh.getParam(prefix + "extrin_calib/Rcl", rotation) ||
                     !nh.getParam(prefix + "extrin_calib/Pcl", translation)))
      throw std::runtime_error("Missing LiDAR2Camera extrinsics: " + id);
    if (rotation.size() != 9 || translation.size() != 3)
      throw std::runtime_error("Invalid camera extrinsic dimensions: " + id);
    view->setImuToLidarExtrinsic(extT, extR);
    view->setLidarToCameraExtrinsic(rotation, translation);
    if (!view->Rcl.allFinite() || !view->Pcl.allFinite() ||
        (view->Rcl.transpose()*view->Rcl-M3D::Identity()).norm() > 1e-4 || view->Rcl.determinant() < .99)
      throw std::runtime_error("Invalid camera extrinsic rotation: " + id);
    view->state = &_state; view->state_propagat = &state_propagat;
    // Every view uses the same kernels with its own configured parameters.
    nh.param<int>(prefix + "grid_size", view->grid_size, grid_size);
    nh.param<int>(prefix + "grid_n_height", view->grid_n_height, grid_n_height);
    nh.param<int>(prefix + "patch_size", view->patch_size, patch_size);
    nh.param<int>(prefix + "patch_pyrimid_level", view->patch_pyrimid_level, patch_pyrimid_level);
    nh.param<int>(prefix + "max_iterations", view->max_iterations, max_iterations);
    nh.param<int>(prefix + "max_observations", view->max_observations, 0);
    nh.param<double>(prefix + "img_point_cov", view->img_point_cov, IMG_POINT_COV);
    nh.param<double>(prefix + "outlier_threshold", view->outlier_threshold, outlier_threshold);
    nh.param<double>(prefix + "time_offset_s", view->time_offset, supplied ? 0 : img_time_offset);
    nh.param<double>(prefix + "exposure_time_offset_s", view->exposure_time_offset, supplied ? 0 : exposure_time_init);
    view->normal_en = normal_en; view->inverse_composition_en = inverse_composition_en;
    view->raycast_en = raycast_en; view->exposure_estimate_en = exposure_estimate_en;
    view->colmap_output_en = colmap_output_en;
    view->plot_flag = false;
    if (view->patch_size <= 0 || view->patch_size % 2 || view->patch_pyrimid_level <= 0 || view->patch_pyrimid_level > 8 ||
        view->grid_size <= 0 || view->grid_n_height <= 0 || view->max_iterations <= 0 || view->max_observations < 0 ||
        !std::isfinite(view->img_point_cov) || view->img_point_cov <= 0 || !std::isfinite(view->time_offset) ||
        !std::isfinite(view->exposure_time_offset)) throw std::runtime_error("Invalid camera settings: " + id);
    view->initializeVIO();
    if (id == anchor) vio_manager->anchor_index = k;
    ROS_INFO("[camera] %s index=%zu topic=%s", id.c_str(), k, view->topic.c_str());
    vio_manager->cameras.push_back(std::move(view));
  }
  if (!ids.count(anchor)) throw std::runtime_error("anchor_camera is not in camera_order");
  for (auto *state : {&_state, &state_propagat, &imu_propagate, &latest_ekf_state}) state->configureCameras(order.size());
  vio_manager->state = &_state; vio_manager->state_propagat = &state_propagat;
  camera_queues.resize(order.size());
  if (img_save_en)
  {
    camera_frames_log.open(string(ROOT_DIR) + "Log/camera_frames.jsonl");
    if (!camera_frames_log) throw std::runtime_error("Cannot open camera frame export");
    for (const auto &camera : vio_manager->cameras)
      std::filesystem::create_directories(string(ROOT_DIR) + "Log/image/" + camera->id);
  }
  vio_manager->audit.open(string(ROOT_DIR) + "Log/camera_observations.csv");
  if (!vio_manager->audit) throw std::runtime_error("Cannot open camera observation audit");
  vio_manager->audit << "timestamp,camera,retrieved,used,inv_exposure,covariance_trace\n";
}

void LIVMapper::img_cbk(const sensor_msgs::ImageConstPtr &msg, size_t camera_index)
{
  if (!img_en) return;
  auto &queue = camera_queues.at(camera_index);
  auto &camera = *vio_manager->cameras.at(camera_index);
  if (hilti_en && ++queue.subsample_counter % 4 != 0) return;
  const double time = msg->header.stamp.toSec() + camera.time_offset + camera.exposure_time_offset;
  if (!std::isfinite(time)) throw std::runtime_error("Non-finite image time");
  if (abs(time - queue.last_timestamp) < .001 || last_timestamp_lidar < 0) return;
  if (time < queue.last_timestamp) throw std::runtime_error("Image time went backwards: " + camera.id);
  if (time - queue.last_timestamp < .02) return;
  std::lock_guard<std::mutex> lock(mtx_buffer);
  queue.images.push_back(getImageFromMsg(msg));
  queue.times.push_back(time);
  queue.stamps.push_back(msg->header.stamp.toNSec() + int64_t(std::llround((camera.time_offset + camera.exposure_time_offset) * 1e9)));
  queue.last_timestamp = time;
  sig_buffer.notify_all();
}

bool LIVMapper::prepareImageBundle(double timestamp)
{
  std::lock_guard<std::mutex> lock(mtx_buffer);
  bool finished = false;
  ros::param::get("/fast_livo/input_finished", finished);
  vector<double> first, last;
  for (auto &queue : camera_queues)
  {
    while (!queue.times.empty() && queue.times.front() < timestamp - vio_manager->sync_tolerance)
    { queue.times.pop_front(); queue.images.pop_front(); queue.stamps.pop_front(); }
    first.push_back(queue.times.empty() ? std::numeric_limits<double>::infinity() : queue.times.front());
    last.push_back(queue.last_timestamp);
  }
  if (!cameraBundleReady(timestamp, first, last, vio_manager->sync_tolerance, camera_maximum_lateness, finished))
    return false;
  vector<MeasureGroup::Image> images;
  for (size_t k = 0; k < camera_queues.size(); ++k)
  {
    auto &queue = camera_queues[k];
    if (!queue.times.empty() && abs(queue.times.front() - timestamp) <= vio_manager->sync_tolerance)
      images.push_back({k, queue.times.front(), queue.images.front(), queue.stamps.front()});
  }
  pending_images = std::move(images);
  return true;
}

void LIVMapper::consumeImageBundle()
{
  for (const auto &image : pending_images)
  {
    auto &queue = camera_queues.at(image.camera_index);
    if (queue.times.empty() || queue.times.front() != image.timestamp)
      throw std::runtime_error("Finalized camera bundle changed before consumption");
    queue.times.pop_front(); queue.images.pop_front(); queue.stamps.pop_front();
  }
  pending_images.clear();
}

void LIVMapper::saveCameraImages(const string &stamp, double timestamp)
{
  camera_frames_log << "{\"timestamp_ns\":" << (std::stoll(stamp.substr(0, stamp.find('.'))) * 1000000000LL + std::stoll(stamp.substr(stamp.find('.') + 1)) * 1000)
                    << ",\"images\":[";
  bool first = true;
  for (const auto &camera : vio_manager->cameras)
  {
    if (!camera->active) continue;
    const string relative = "image/" + (camera->camera_index == vio_manager->anchor_index ? "" : camera->id + "/") + stamp + ".png";
    if (!cv::imwrite(string(ROOT_DIR) + "Log/" + relative, camera->img_rgb))
      throw std::runtime_error("Failed to save camera image");
    SE3 T_W_C = camera->new_frame_->T_f_w_.inverse();
    Eigen::Quaterniond q(T_W_C.rotation_matrix());
    if (!first) camera_frames_log << ',';
    first = false;
    camera_frames_log << std::setprecision(17) << "{\"camera_id\":\"" << camera->id
      << "\",\"timestamp_ns\":" << camera->image_timestamp_ns
      << ",\"state_timestamp_ns\":" << int64_t(std::llround(timestamp * 1e9))
      << ",\"path\":\"raw/" << relative << "\",\"T_W_C\":{\"translation_m\":["
      << T_W_C.translation().x() << ',' << T_W_C.translation().y() << ',' << T_W_C.translation().z()
      << "],\"quaternion_xyzw\":[" << q.x() << ',' << q.y() << ',' << q.z() << ',' << q.w() << "]}}";
  }
  camera_frames_log << "]}\n" << std::flush;
}
