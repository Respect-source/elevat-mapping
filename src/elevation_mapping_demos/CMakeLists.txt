cmake_minimum_required(VERSION 3.8)
project(elevation_mapping_demos)

# Default to C++14 (ROS2 推荐)
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

# Find dependencies
find_package(ament_cmake REQUIRED)
find_package(elevation_mapping REQUIRED)
find_package(grid_map_rviz_plugin REQUIRED)
find_package(rviz_common REQUIRED)
find_package(rviz_rendering REQUIRED)

# Install resource files
install(DIRECTORY config doc launch rviz sample_data
  DESTINATION share/${PROJECT_NAME}
)

# Install Python scripts
install(PROGRAMS scripts/tf_to_pose_publisher.py
  DESTINATION lib/${PROJECT_NAME}
)

# Tests
if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  ament_lint_auto_find_test_dependencies()

  find_package(ament_cmake_gtest REQUIRED)
  ament_add_gtest(test_${PROJECT_NAME} test/empty_test.cpp)
  if(TARGET test_${PROJECT_NAME})
    target_include_directories(test_${PROJECT_NAME}
      PRIVATE
        include
    )
    ament_target_dependencies(test_${PROJECT_NAME}
      elevation_mapping
      grid_map_rviz_plugin
      rviz_common
      rviz_rendering
    )
  endif()
endif()

ament_package()
