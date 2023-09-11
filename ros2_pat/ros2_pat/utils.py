import rclpy
import rclpy.node, rclpy.clock
import tf_transformations
import tf2_ros
from tf2_ros.static_transform_broadcaster import StaticTransformBroadcaster
from geometry_msgs.msg import TransformStamped

import numpy as np

class PatTransformPublisher:
    def __init__(self, parent_node):
        self.tf_broadcaster = tf2_ros.TransformBroadcaster(parent_node)
        self.tf_static_broadcaster = StaticTransformBroadcaster(parent_node)
        #self.send_w_enu_to_ned_transform(rclpy.clock.Clock().now())

    def send_w_enu_to_ned_transform(self, t):
        R_enu2ned = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        T_enu2ned = np.eye(4); T_enu2ned[:3,:3] = R_enu2ned
        self.send_transform("w_enu", "w_ned", t, T_enu2ned, static=True)

    def send_w_ned_to_b_transform(self, t, T_w2b):
        self.send_transform("w_ned", "b_frd", t, T_w2b)
        
    def send_transform(self, f1, f2, t, T_f1tof2, static=False):
        tf_msg = TransformStamped()
        tf_msg.header.frame_id = f1
        tf_msg.child_frame_id = f2
        tf_msg.header.stamp = t.to_msg()
        _r = tf_msg.transform.rotation
        _r.x, _r.y, _r.z, _r.w = tf_transformations.quaternion_from_matrix(T_f1tof2)
        _t = tf_msg.transform.translation
        _t.x, _t.y, _t.z = T_f1tof2[:3,3]
        if static:
            self.tf_static_broadcaster.sendTransform(tf_msg)
        else:
            self.tf_broadcaster.sendTransform(tf_msg)


        
