#!/bin/bash

echo "=== 高程图系统诊断脚本 ==="
echo ""

# 设置环境
cd /home/panqiulong/rc_2026_nav
source install/setup.bash

echo "1. 检查ROS域ID:"
echo "ROS_DOMAIN_ID = $ROS_DOMAIN_ID"
echo ""

echo "2. 检查运行中的节点:"
ros2 node list
echo ""

echo "3. 检查话题列表:"
ros2 topic list
echo ""

echo "4. 检查点云数据:"
echo "检查 /livox/lidar 话题:"
ros2 topic info /livox/lidar
echo ""

echo "5. 检查IMU数据:"
echo "检查 /livox/imu 话题:"
ros2 topic info /livox/imu
echo ""

echo "6. 检查TF树:"
ros2 run tf2_tools view_frames
echo ""

echo "7. 检查TF变换:"
echo "检查 camera_init -> aft_mapped:"
timeout 5 ros2 run tf2_ros tf2_echo camera_init aft_mapped || echo "TF变换不可用"
echo ""

echo "检查 aft_mapped -> livox_frame:"
timeout 5 ros2 run tf2_ros tf2_echo aft_mapped livox_frame || echo "TF变换不可用"
echo ""

echo "8. 检查高程图话题:"
echo "检查 /elevation_map 话题:"
ros2 topic info /elevation_map
echo ""

echo "9. 检查高程图节点参数:"
ros2 param list /elevation_mapping 2>/dev/null || echo "高程图节点未运行"
echo ""

echo "10. 检查点云数据内容:"
echo "获取一次点云数据:"
timeout 5 ros2 topic echo /livox/lidar --once || echo "无法获取点云数据"
echo ""

echo "诊断完成！"
