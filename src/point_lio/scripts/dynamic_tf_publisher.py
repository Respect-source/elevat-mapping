#!/usr/bin/env python3
import rclpy
from rclpy.node import Node
from geometry_msgs.msg import TransformStamped
from tf2_ros import TransformBroadcaster

class DynamicTFPublisher(Node):
    def __init__(self):
        super().__init__('dynamic_tf_publisher')
        self.tf_broadcaster = TransformBroadcaster(self)

        # 固定外参：aft_mapped -> livox_frame（请按实际安装修改）
        self.extrinsic_tx = 0.15
        self.extrinsic_ty = 0.0
        self.extrinsic_tz = 0.3

        # 不再发布 camera_init->aft_mapped 或 aft_mapped->body 的静态 TF，避免与里程计发布冲突

        # 动态定时发布：aft_mapped -> livox_frame（仅位置外参，无姿态，姿态由上游 odom 提供到 aft_mapped）
        self.timer = self.create_timer(0.05, self.publish_lidar_extrinsic)
        self.get_logger().info('TF 发布器已启动：仅发布 aft_mapped->livox_frame 外参')

    def publish_lidar_extrinsic(self):
        stamp = self.get_clock().now().to_msg()
        t = TransformStamped()
        t.header.stamp = stamp
        t.header.frame_id = 'aft_mapped'
        t.child_frame_id = 'livox_frame'
        t.transform.translation.x = self.extrinsic_tx
        t.transform.translation.y = self.extrinsic_ty
        t.transform.translation.z = self.extrinsic_tz
        t.transform.rotation.x = 0.0
        t.transform.rotation.y = 0.0
        t.transform.rotation.z = 0.0
        t.transform.rotation.w = 1.0
        self.tf_broadcaster.sendTransform(t)


def main(args=None):
    rclpy.init(args=args)
    node = DynamicTFPublisher()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        node.get_logger().info('节点已停止')
    finally:
        node.destroy_node()
        rclpy.shutdown()

if __name__ == '__main__':
    main()