##
#
# python3  ros2_pat/ros2_pat/traj_publisher.py --ros-args -p traj_name:=si2 -p list:=true



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

FOO = True

class TrajPublisher(Node):

    def __init__(self):
        super().__init__('pat3_trajectory_publisher')

        self.declare_parameter('traj_name', 'si1')
        self.declare_parameter('list', False)

        if (self.get_parameter('list').get_parameter_value().bool_value):
             all_traj = '\n'.join(trjf.list())
             self.get_logger().info(f"available trajectories: \n{all_traj}")
        
        traj_name = self.get_parameter('traj_name').get_parameter_value().string_value
        self.get_logger().info(f"loading trajectory: {traj_name}")
        #self.traj = trjf.Traj43(duration=7.)
        self.traj, desc = trjf.get(traj_name)
        self.get_logger().info(f'  Description: {desc}')
        
        self.tf_pub = ros2_pat.utils.TransformPublisher(self)
        self.tf_pub.send_w_enu_to_ned_transform(rclpy.clock.Clock().now())

        

        if FOO:
            self.marker_pub = self.create_publisher(Marker, "/pat/reference", 2)
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
            #self.marker.mesh_resource = "package://ros2_pat/meshes/quadcopter.dae"
            self.marker.mesh_resource = "package://ros2_pat/meshes/quad.dae"
            #self.marker.mesh_resource = "package://ros2_pat/meshes/mavic.dae"
        else:
            self.marker_pub = ros2_pat.utils.MarkerPublisher(self, '/pat/reference', mtype=10,
                                                             mesh='package://ros2_pat/meshes/quad.dae')


        try:
            nwp = len(self.traj._g.waypoints)
            self.wp_pub = ros2_pat.utils.MarkerArrayPublisher(self, 'pat/waypoints', mtype=2,
                                                              meshes=['']*nwp,
                                                              colors=[[1., 0.2, 0.2, 0.5]]*nwp, scales=[(0.1, 0.1, 0.1)]*nwp)
        except AttributeError: pass
            
        self.trj_pub = ros2_pat.utils.TrajPublisher(self.traj, self, 'pat/ref_trajectory')

            
        timer_period = 1./20.  # seconds
        self.timer = self.create_timer(timer_period, self.timer_callback)

        
        self.df = ctl.DiffFlatness()
        self.fdm = fdm.FDM()
        self.t0 = Clock().now().nanoseconds*1e-9

        

    def timer_callback(self):
        now = Clock().now()
        t = np.fmod(now.nanoseconds*1e-9 - self.t0, self.traj.duration)
        Y = self.traj.get(t)
        X, U, Xd = self.df.state_and_cmd_of_flat_output(Y, self.fdm.P)
        if FOO:
            self.marker.header.stamp = now.to_msg()
            self.marker.pose.position.x = X[fdm.sv_x] #np.sin(self.i%100/50.*np.pi) #0.
            self.marker.pose.position.y = X[fdm.sv_y]
            self.marker.pose.position.z = X[fdm.sv_z]
            self.marker.pose.orientation.x = X[fdm.sv_qx]
            self.marker.pose.orientation.y = X[fdm.sv_qy]
            self.marker.pose.orientation.z = X[fdm.sv_qz]
            self.marker.pose.orientation.w = X[fdm.sv_qi]
            self.marker_pub.publish(self.marker)
        else:
            #T_ned2b = np.eye(4)
            T_ned2b = ros2_pat.utils.T_of_t_q(X[fdm.sv_slice_pos], X[fdm.sv_slice_quat])
            self.marker_pub.publish(now, T_ned2b)
        
        self.trj_pub.publish()
        try:
            nwp = len(self.traj._g.waypoints)
            T_ned2b = np.array([np.eye(4)]*nwp)
            for i in range(nwp):
                T_ned2b[i][:3,3] = self.traj._g.waypoints[i]
            self.wp_pub.publish(T_ned2b)
        except AttributeError: pass 
        
        
def main(args=None):
    rclpy.init(args=args)
    traj_publisher = TrajPublisher()
    rclpy.spin(traj_publisher)
    traj_publisher.destroy_node()
    rclpy.shutdown()


if __name__ == '__main__':
    main()
