#!/usr/bin/env python3
import rclpy
from rclpy.node import Node
from geometry_msgs.msg import TransformStamped
from tf2_ros import TransformBroadcaster
from sensor_msgs.msg import Imu

class DynamicTFPublisher(Node):
    def __init__(self):
        super().__init__('dynamic_tf_publisher')
        self.tf_broadcaster = TransformBroadcaster(self)
        self.imu_sub = self.create_subscription(Imu, '/livox/imu', self.imu_callback, 10)
        self.get_logger().info("动态 TF 发布节点已启动，订阅 /livox/imu")

    def imu_callback(self, msg):
        # 发布 aft_mapped -> livox_frame
        t = TransformStamped()
        t.header.stamp = msg.header.stamp
        t.header.frame_id = 'aft_mapped'
        t.child_frame_id = 'livox_frame'
        t.transform.translation.x = 0.0  # 根据激光雷达安装位置调整
        t.transform.translation.y = 0.0
        t.transform.translation.z = 0.0
        t.transform.rotation = msg.orientation  # 使用 IMU 四元数
        self.tf_broadcaster.sendTransform(t)
        self.get_logger().info(f"发布 TF: aft_mapped -> livox_frame, stamp={t.header.stamp}")

        # 临时发布 camera_init -> aft_mapped（占位符）
        t2 = TransformStamped()
        t2.header.stamp = msg.header.stamp
        t2.header.frame_id = 'camera_init'
        t2.child_frame_id = 'aft_mapped'
        t2.transform.translation.x = 0.0
        t2.transform.translation.y = 0.0
        t2.transform.translation.z = 0.0
        t2.transform.rotation.x = 0.0
        t2.transform.rotation.y = 0.0
        t2.transform.rotation.z = 0.0
        t2.transform.rotation.w = 1.0
        self.tf_broadcaster.sendTransform(t2)
        self.get_logger().info(f"发布 TF: camera_init -> aft_mapped, stamp={t2.header.stamp}")

def main(args=None):
    rclpy.init(args=args)
    node = DynamicTFPublisher()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        node.get_logger().info("节点已停止")
    finally:
        node.destroy_node()
        rclpy.shutdown()

if __name__ == '__main__':
    main()