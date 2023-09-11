import rclpy
from rclpy.node import Node
from rclpy.clock import Clock

from std_msgs.msg import String
from visualization_msgs.msg import Marker
import geometry_msgs.msg

import tf_transformations
import tf2_ros

import numpy as np, math

import pat3.vehicles.rotorcraft.multirotor_trajectory as trj
import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as trjf
import pat3.vehicles.rotorcraft.multirotor_control as ctl
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.algebra as pal

import ros2_pat.utils

class TrajPublisher(Node):

    def __init__(self):
        super().__init__('minimal_publisher')

        self.tf_pub = ros2_pat.utils.PatTransformPublisher(self)
        self.tf_pub.send_w_enu_to_ned_transform(rclpy.clock.Clock().now())


        self.marker_pub = self.create_publisher(Marker, "/visualization_marker", 2)
        self.marker = Marker()
        self.marker.header.frame_id = "/w_ned"
        self.marker.type = 2
        self.marker.id = 0
        self.marker.scale.x = 1.0
        self.marker.scale.y = 1.0
        self.marker.scale.z = 1.0
        self.marker.color.r = 0.0
        self.marker.color.g = 1.0
        self.marker.color.b = 0.0
        self.marker.color.a = 1.0
        self.marker.type = 10 #visualization_msgs::Marker::MESH_RESOURCE;
        self.marker.mesh_resource = "package://ros2_pat/meshes/quadcopter.dae"
        #self.marker.mesh_resource = "package://ros2_pat/meshes/quad.dae"
        #self.marker.mesh_resource = "package://ros2_pat/meshes/mavic.dae"

        timer_period = 1./20.  # seconds
        self.timer = self.create_timer(timer_period, self.timer_callback)
        #self.i = 0
        #self.traj = trj.Circle(v=2.)
        #self.traj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[1, 0, 0, 0], duration=2.)
        #self.traj = trj.SmoothBackAndForth(x0=[0, 0, 0, np.pi/2], x1=[1, 0, 0, np.pi/2])
        #self.traj = trj.SmoothBackAndForth(x0=[0, 0, 0, 0], x1=[1, 0, 0, 0])
        #self.traj = trjf.Traj14()
        self.traj = trjf.Traj43()
        self.df = ctl.DiffFlatness()
        self.fdm = fdm.FDM()
        self.t0 = Clock().now().nanoseconds*1e-9

        

    def timer_callback(self):
        now = Clock().now()
        t = np.fmod(now.nanoseconds*1e-9 - self.t0, self.traj.duration)
        Y = self.traj.get(t)
        X, U, Xd = self.df.state_and_cmd_of_flat_output(Y, self.fdm.P)
        
        self.marker.header.stamp = now.to_msg()
        self.marker.pose.position.x = X[fdm.sv_x] #np.sin(self.i%100/50.*np.pi) #0.
        self.marker.pose.position.y = X[fdm.sv_y]
        self.marker.pose.position.z = X[fdm.sv_z]
        self.marker.pose.orientation.x = X[fdm.sv_qx]
        self.marker.pose.orientation.y = X[fdm.sv_qy]
        self.marker.pose.orientation.z = X[fdm.sv_qz]
        self.marker.pose.orientation.w = X[fdm.sv_qi]
        self.marker_pub.publish(self.marker)



        
        
def main(args=None):
    rclpy.init(args=args)
    traj_publisher = TrajPublisher()
    rclpy.spin(traj_publisher)
    traj_publisher.destroy_node()
    rclpy.shutdown()


if __name__ == '__main__':
    main()
