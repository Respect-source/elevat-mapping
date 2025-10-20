#!/usr/bin/env python3

import rclpy
from rclpy.node import Node
from grid_map_msgs.msg import GridMap
from sensor_msgs.msg import PointCloud2
import time

class ElevationMappingTester(Node):
    def __init__(self):
        super().__init__('elevation_mapping_tester')
        
        # 订阅高程图
        self.elevation_map_sub = self.create_subscription(
            GridMap,
            '/elevation_map',
            self.elevation_map_callback,
            10
        )
        
        # 订阅点云数据
        self.pointcloud_sub = self.create_subscription(
            PointCloud2,
            '/livox/lidar',
            self.pointcloud_callback,
            10
        )
        
        self.elevation_map_received = False
        self.pointcloud_received = False
        self.elevation_map_count = 0
        self.pointcloud_count = 0
        
        self.get_logger().info('高程图测试节点已启动')
        self.get_logger().info('等待数据...')

    def elevation_map_callback(self, msg):
        self.elevation_map_received = True
        self.elevation_map_count += 1
        
        self.get_logger().info(f'收到高程图 #{self.elevation_map_count}:')
        self.get_logger().info(f'  时间戳: {msg.header.stamp.sec}.{msg.header.stamp.nanosec}')
        self.get_logger().info(f'  坐标系: {msg.header.frame_id}')
        self.get_logger().info(f'  层数: {len(msg.layers)}')
        
        for layer in msg.layers:
            self.get_logger().info(f'    层: {layer}')
        
        self.get_logger().info(f'  地图大小: {msg.info.length_x:.2f} x {msg.info.length_y:.2f}')
        self.get_logger().info(f'  分辨率: {msg.info.resolution:.3f}')
        self.get_logger().info('---')

    def pointcloud_callback(self, msg):
        if not self.pointcloud_received:
            self.pointcloud_received = True
            self.get_logger().info('收到点云数据!')
        
        self.pointcloud_count += 1
        if self.pointcloud_count % 10 == 0:  # 每10个点云打印一次
            self.get_logger().info(f'已收到 {self.pointcloud_count} 个点云数据')

    def check_status(self):
        if self.pointcloud_received and not self.elevation_map_received:
            self.get_logger().warn('点云数据正常，但高程图没有输出！')
        elif not self.pointcloud_received:
            self.get_logger().warn('没有收到点云数据！')
        elif self.elevation_map_received:
            self.get_logger().info('✅ 高程图系统工作正常！')

def main(args=None):
    rclpy.init(args=args)
    node = ElevationMappingTester()
    
    try:
        # 运行30秒
        start_time = time.time()
        while time.time() - start_time < 30:
            rclpy.spin_once(node, timeout_sec=1.0)
            if time.time() - start_time > 10:  # 10秒后检查状态
                node.check_status()
                break
    except KeyboardInterrupt:
        node.get_logger().info('测试被中断')
    finally:
        node.destroy_node()
        rclpy.shutdown()

if __name__ == '__main__':
    main()
