import rclpy
import rclpy.node, rclpy.clock
import tf_transformations

import tf2_ros
from tf2_ros.static_transform_broadcaster import StaticTransformBroadcaster
from geometry_msgs.msg import TransformStamped

import numpy as np

import ros2_pat.utils

class TestNode(rclpy.node.Node):
    def __init__(self):
        super().__init__('test_orientation_node')
        self.tf_pub = ros2_pat.utils.PatTransformPublisher(self)
        self.tf_pub.send_w_enu_to_ned_transform(rclpy.clock.Clock().now())

        timer_period = 1./20.  # seconds
        self.timer = self.create_timer(timer_period, self.timer_callback)

    def timer_callback(self):
        phi, theta, psi = 0, 0, np.deg2rad(10.)
        T_ned_to_body = tf_transformations.euler_matrix(phi, theta, psi, axes='sxyz')
        T_ned_to_body[:3,3] = [1,0,0]
        #print(T_ned_to_body)
        self.tf_pub.send_w_ned_to_b_transform(rclpy.clock.Clock().now(), T_ned_to_body)
        
def main(args=None):
    rclpy.init(args=args)
    test_node = TestNode()
    rclpy.spin(test_node)
    rclpy.shutdown()


if __name__ == '__main__':
    main()
