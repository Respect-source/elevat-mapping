import os
from launch import LaunchDescription
from launch_ros.actions import Node
from ament_index_python.packages import get_package_share_directory
from launch.actions import IncludeLaunchDescription
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import PathJoinSubstitution


def generate_launch_description():
    config = os.path.join(
        get_package_share_directory("elevation_mapping"),
        "config",
        "sensor_processors",
        "mapping.yaml"
    )

    rviz_config = os.path.join(
        get_package_share_directory("elevation_mapping"),
        "rviz",
        "elevation_mapping.rviz"   # 这里假设你有 rviz 配置文件
    )

    point_lio_launch = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            PathJoinSubstitution([
                get_package_share_directory("point_lio"),
                "launch",
                "point_lio.launch.py",
            ])
        ),
        launch_arguments={
            "rviz": "False",
            # 如需自定义 point_lio 配置，取消注释并填写路径
            # "point_lio_cfg_dir": "/absolute/path/to/your/mid360.yaml",
        }.items(),
    )

    return LaunchDescription([
        #point_lio_launch,
        Node(
            package="elevation_mapping",
            executable="elevation_mapping_node",
            name="elevation_mapping",
            output="screen",
            parameters=[config],
        ),
    ])
