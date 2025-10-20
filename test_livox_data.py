#!/usr/bin/env python3

import rclpy
from rclpy.node import Node
from sensor_msgs.msg import PointCloud2
import struct

class LivoxDataChecker(Node):
    def __init__(self):
        super().__init__('livox_data_checker')
        self.subscription = self.create_subscription(
            PointCloud2,
            '/livox/lidar',
            self.pointcloud_callback,
            10
        )
        self.get_logger().info('等待点云数据...')

    def pointcloud_callback(self, msg):
        self.get_logger().info(f'收到点云数据:')
        self.get_logger().info(f'  时间戳: {msg.header.stamp.sec}.{msg.header.stamp.nanosec}')
        self.get_logger().info(f'  坐标系: {msg.header.frame_id}')
        self.get_logger().info(f'  点云大小: {msg.width} x {msg.height}')
        self.get_logger().info(f'  字段数量: {len(msg.fields)}')
        
        # 打印字段信息
        for i, field in enumerate(msg.fields):
            self.get_logger().info(f'  字段 {i}: {field.name} (类型: {field.datatype}, 偏移: {field.offset})')
        
        # 检查是否有时间戳字段
        has_timestamp = any(field.name in ['t', 'time', 'timestamp'] for field in msg.fields)
        self.get_logger().info(f'  包含时间戳字段: {has_timestamp}')
        
        # 检查点云数据的前几个点
        if msg.data:
            point_size = msg.point_step
            self.get_logger().info(f'  点大小: {point_size} 字节')
            
            # 解析前3个点
            for i in range(min(3, msg.width)):
                offset = i * point_size
                point_data = msg.data[offset:offset + point_size]
                
                # 解析 x, y, z (假设是前3个float32)
                if len(point_data) >= 12:
                    x, y, z = struct.unpack('fff', point_data[:12])
                    self.get_logger().info(f'  点 {i}: x={x:.3f}, y={y:.3f}, z={z:.3f}')
                
                # 检查时间戳字段
                if has_timestamp:
                    for field in msg.fields:
                        if field.name in ['t', 'time', 'timestamp']:
                            if field.datatype == 7:  # UINT64
                                t_offset = field.offset
                                if len(point_data) > t_offset + 8:
                                    timestamp = struct.unpack('Q', point_data[t_offset:t_offset + 8])[0]
                                    self.get_logger().info(f'  时间戳: {timestamp} ns')
                            break
        
        self.get_logger().info('---')

def main(args=None):
    rclpy.init(args=args)
    node = LivoxDataChecker()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        node.get_logger().info('节点已停止')
    finally:
        node.destroy_node()
        rclpy.shutdown()

if __name__ == '__main__':
    main()
