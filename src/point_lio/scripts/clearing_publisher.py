#!/usr/bin/env python3
import rclpy
from rclpy.node import Node
from sensor_msgs.msg import PointCloud2, PointField
from std_msgs.msg import Header
import struct

class ClearingPublisher(Node):
    def __init__(self):
        super().__init__('clearing_publisher')
        self.publisher_ = self.create_publisher(PointCloud2, '/elevation_map_clearing', 10)
        self.timer = self.create_timer(0.1, self.publish_clearing)  # 10Hz

    def publish_clearing(self):
        msg = PointCloud2()
        msg.header = Header()
        msg.header.frame_id = 'aft_mapped'  # 与mapping.yaml的map_frame_id一致
        msg.header.stamp = self.get_clock().now().to_msg()
        resolution = 0.05  # 匹配mapping.yaml的resolution
        size = 12.0  # 匹配length_in_x, length_in_y
        grid_size = int(size / resolution)  # 240x240
        msg.height = 1
        msg.width = grid_size * grid_size  # 57,600点
        msg.fields = [
            PointField(name='x', offset=0, datatype=PointField.FLOAT32, count=1),
            PointField(name='y', offset=4, datatype=PointField.FLOAT32, count=1),
            PointField(name='z', offset=8, datatype=PointField.FLOAT32, count=1)
        ]
        msg.is_bigendian = False
        msg.point_step = 12
        msg.row_step = msg.point_step * msg.width
        msg.is_dense = False
        data = []
        for i in range(grid_size):
            for j in range(grid_size):
                x = -size / 2 + i * resolution
                y = -size / 2 + j * resolution
                data.extend(struct.pack('<fff', x, y, float('nan')))
        msg.data = bytes(data)
        self.publisher_.publish(msg)
        self.get_logger().info(f'Published clearing message with {msg.width} NaN points')

def main(args=None):
    rclpy.init(args=args)
    node = ClearingPublisher()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()