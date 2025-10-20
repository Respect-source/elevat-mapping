import os
from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.conditions import IfCondition
from launch.substitutions import LaunchConfiguration, PathJoinSubstitution
from launch_ros.actions import Node
from launch.actions import IncludeLaunchDescription
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import PathJoinSubstitution
from launch.actions import SetEnvironmentVariable



def generate_launch_description():
    # Map fully qualified names to relative ones so the node's namespace can be prepended.
    # In case of the transforms (tf), currently, there doesn't seem to be a better alternative
    # https://github.com/ros/geometry2/issues/32
    # https://github.com/ros/robot_state_publisher/pull/30
    # TODO(orduno) Substitute with `PushNodeRemapping`
    #              https://github.com/ros2/launch_ros/issues/56
    
    remappings = [("/tf", "tf"), ("/tf_static", "tf_static")]

    namespace = LaunchConfiguration("namespace")
    use_rviz = LaunchConfiguration("rviz")
    point_lio_cfg_dir = LaunchConfiguration("point_lio_cfg_dir")

    point_lio_dir = get_package_share_directory("point_lio")

    declare_namespace = DeclareLaunchArgument(
        "namespace",
        default_value="",
        description="Namespace for the node",
    )

    declare_rviz = DeclareLaunchArgument(
        "rviz", default_value="True", description="Flag to launch RViz."
    )

    declare_point_lio_cfg_dir = DeclareLaunchArgument(
        "point_lio_cfg_dir",
        default_value=PathJoinSubstitution([point_lio_dir, "config", "mid360.yaml"]),
        description="Path to the Point-LIO config file",
    )

    start_point_lio_node = Node(
        package="point_lio",
        executable="pointlio_mapping",
        namespace=namespace,
        parameters=[point_lio_cfg_dir],
        remappings=remappings,
        output="screen",
    )

    start_rviz_node = Node(
        condition=IfCondition(use_rviz),
        package="rviz2",
        executable="rviz2",
        namespace=namespace,
        name="rviz",
        remappings=remappings,
        arguments=[
            "-d",
            PathJoinSubstitution([point_lio_dir, "rviz_cfg", "loam_livox"]),
            ".rviz",
        ],
    )
    
    elevation_mapping = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            PathJoinSubstitution([
                get_package_share_directory("elevation_mapping"),
                "launch",
                "elevation_mapping.launch.py",
            ])
        )
    )

    msg_laser = IncludeLaunchDescription(
    PythonLaunchDescriptionSource(
        PathJoinSubstitution([
            get_package_share_directory("livox_ros_driver2"),
            "launch_ROS2",
            "msg_MID360_launch.py",
        ])
    ),
    launch_arguments={'publish_custom_msg': 'false'}.items()
)
    
    # 静态 TF：camera_init -> aft_mapped
    #static_tf_1 = Node(
    #   package="tf2_ros", executable="static_transform_publisher",
    #    arguments=["--x","0","--y","0","--z","0","--qx","0","--qy","0","--qz","0","--qw","1",
    #           "--frame-id","camera_init","--child-frame-id","aft_mapped"],
    #    output="screen",
    #)

    # 静态 TF：aft_mapped -> livox_frame
    #static_tf_2 = Node(
    #    package="tf2_ros", executable="static_transform_publisher",
    #    arguments=["--x","0","--y","0","--z","0","--qx","0","--qy","0","--qz","0","--qw","1",
    #           "--frame-id","aft_mapped","--child-frame-id","livox_frame"],
    #    output="screen",
    #)
    dynamic_tf_node = Node(
       package="point_lio",  # 如果你的包名不同，请替换
       executable="dynamic_tf_publisher.py",
       name="dynamic_tf_publisher",
       output="screen",
    )
    #clearing_publisher_node = Node(
    #package="point_lio",
    #executable="clearing_publisher.py",
    #name="clearing_publisher",
    #output="screen",
#)
    ld = LaunchDescription()

    ld.add_action(SetEnvironmentVariable('RMW_IMPLEMENTATION', 'rmw_cyclonedds_cpp'))
    ld.add_action(declare_namespace)
    ld.add_action(declare_rviz)
    ld.add_action(msg_laser)
    ld.add_action(declare_point_lio_cfg_dir)
    ld.add_action(start_point_lio_node)
    ld.add_action(start_rviz_node)
    ld.add_action(elevation_mapping)
    #ld.add_action(static_tf_1)
    #ld.add_action(static_tf_2)
    ld.add_action(dynamic_tf_node)
    #ld.add_action(clearing_publisher_node)
    return ld
