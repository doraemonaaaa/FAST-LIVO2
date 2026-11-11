cmake_minimum_required(VERSION 2.8.3)
project(fast_livo)

set(CMAKE_BUILD_TYPE "Release")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set common compile options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -fexceptions")

# Specific settings for Debug build
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")

# Detect CPU architecture
message(STATUS "Current CPU architecture: ${CMAKE_SYSTEM_PROCESSOR}")

# Specific settings for Release build
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm|aarch64|ARM|AARCH64)")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
    # 64-bit ARM optimizations (e.g., RK3588 and Jetson Orin NX)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mcpu=native -mtune=native -ffast-math")
    message(STATUS "Using 64-bit ARM optimizations: -O3 -mcpu=native -mtune=native -ffast-math")
  else()
    # 32-bit ARM optimizations with NEON support
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mcpu=native -mtune=native -mfpu=neon -ffast-math")
    message(STATUS "Using 32-bit ARM optimizations: -O3 -mcpu=native -mtune=native -mfpu=neon -ffast-math")
  endif()
  add_definitions(-DARM_ARCH)
else()
  # x86-64 (Intel/AMD) optimizations
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -mtune=native -funroll-loops") #-flto
  message(STATUS "Using general x86 optimizations: -O3 -march=native -mtune=native -funroll-loops") 
  add_definitions(-DX86_ARCH)
endif()

# Define project root directory
add_definitions(-DROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\")

# Detect CPU core count for potential multithreading optimization
include(ProcessorCount)
ProcessorCount(N)
message(STATUS "Processor count: ${N}")

# Set the number of cores for multithreading
if(N GREATER 4)
  math(EXPR PROC_NUM "4")
  add_definitions(-DMP_EN -DMP_PROC_NUM=${PROC_NUM})
  message(STATUS "Multithreading enabled. Cores: ${PROC_NUM}")
elseif(N GREATER 1)
  math(EXPR PROC_NUM "${N}")
  add_definitions(-DMP_EN -DMP_PROC_NUM=${PROC_NUM})
  message(STATUS "Multithreading enabled. Cores: ${PROC_NUM}")
else()
  add_definitions(-DMP_PROC_NUM=1)
  message(STATUS "Single core detected. Multithreading disabled.")
endif()

# Check for OpenMP support
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found")
  add_compile_options(${OpenMP_CXX_FLAGS})
else()
  message(STATUS "OpenMP not found, proceeding without it")
endif()

# Check for mimalloc support
find_package(mimalloc QUIET)
if(mimalloc_FOUND)
  message(STATUS "mimalloc found")
else()
  message(STATUS "mimalloc not found, proceeding without it")
endif()



# Find catkin and required dependencies
find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  nav_msgs
  sensor_msgs
  roscpp
  rospy
  std_msgs
  pcl_ros
  tf
  message_generation
  eigen_conversions
  vikit_common
  vikit_ros
  cv_bridge
  image_transport
  rosbag
)

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Sophus REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(Ceres REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP REQUIRED jsoncpp)

set(Sophus_LIBRARIES libSophus.so)

# Define the catkin package
catkin_package(
  CATKIN_DEPENDS geometry_msgs nav_msgs roscpp rospy std_msgs message_runtime cv_bridge vikit_common vikit_ros image_transport
  DEPENDS EIGEN3 PCL OpenCV Sophus
)

# Include directories for dependencies
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${PCL_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${Sophus_INCLUDE_DIRS}
  include
  ${CERES_INCLUDE_DIRS}
  ${JSONCPP_INCLUDE_DIRS}
)

# Add libraries
add_library(vio src/vio.cpp src/multi_camera_vio.cpp src/frame.cpp src/visual_point.cpp)
add_library(lio src/voxel_map.cpp)
add_library(pre src/preprocess.cpp)
add_library(initialization
  src/initialization/batch_imu_initializer.cpp
  src/initialization/batch_lidar_trajectory.cpp
  src/initialization/imu_time_offset/lidar_rotation.cpp
  src/initialization/imu_time_offset/scan_matching.cpp
  src/initialization/imu_time_offset/imu_time_offset.cpp)
target_link_libraries(initialization ${CERES_LIBRARIES} ${JSONCPP_LIBRARIES} ${PCL_LIBRARIES} ${catkin_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
add_executable(fastlivo_batch_initialize src/initialization/batch_main.cpp)
target_link_libraries(fastlivo_batch_initialize initialization pre yaml-cpp ${catkin_LIBRARIES} ${PCL_LIBRARIES})
add_library(imu_proc src/IMU_Processing.cpp src/initialization/imu_rewind.cpp)
add_library(laser_mapping src/LIVMapper.cpp src/multi_camera_input.cpp src/sequential_vio.cpp src/scan_vio.cpp src/initialization/startup_replay.cpp src/initialization/staged_initialization.cpp)
target_sources(laser_mapping PRIVATE src/initialization/gravity_norm_constraint.cpp)
target_link_libraries(laser_mapping initialization)

# Add the main executable
add_executable(fastlivo_mapping src/main.cpp)

# Link libraries to the executable
target_link_libraries(fastlivo_mapping
  laser_mapping
  vio
  lio
  pre
  imu_proc
  ${catkin_LIBRARIES}
  ${PCL_LIBRARIES}
  ${OpenCV_LIBRARIES}
  ${Sophus_LIBRARIES}
  ${Boost_LIBRARIES}
)

# Link mimalloc if found
if(mimalloc_FOUND)
  target_link_libraries(fastlivo_mapping mimalloc)
endif()

if(CATKIN_ENABLE_TESTING)
  add_executable(batch_imu_initialization_test test/batch_imu_initialization_test.cpp)
  target_link_libraries(batch_imu_initialization_test initialization ${catkin_LIBRARIES})
  add_test(NAME batch_imu_initialization COMMAND batch_imu_initialization_test)
  add_executable(gravity_norm_constraint_test test/gravity_norm_constraint_test.cpp)
  target_link_libraries(gravity_norm_constraint_test ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
  add_test(NAME gravity_norm_constraint COMMAND gravity_norm_constraint_test)
  add_executable(multi_camera_test test/multi_camera_test.cpp)
  target_link_libraries(multi_camera_test vio lio ${catkin_LIBRARIES} ${PCL_LIBRARIES}
                        ${OpenCV_LIBRARIES} ${Sophus_LIBRARIES} ${Boost_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
  add_test(NAME multi_camera_algebra_and_jacobians COMMAND multi_camera_test)
endif()
