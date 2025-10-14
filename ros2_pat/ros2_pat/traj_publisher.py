##
#
# python3  ros2_pat/ros2_pat/traj_publisher.py --ros-args -p traj_name:=si2 -p list:=true
#


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
        super().__init__('pat3_trajectory_publisher')

        self.declare_parameter('traj_name', 'si1')
        self.declare_parameter('list', False)

        if (self.get_parameter('list').get_parameter_value().bool_value):
             all_traj = '\n'.join(trjf.list())
             self.get_logger().info(f"available trajectories: \n{all_traj}")
        
        traj_name = self.get_parameter('traj_name').get_parameter_value().string_value
        self.get_logger().info(f"loading trajectory: {traj_name}")
        self.traj, desc = trjf.get(traj_name)
        self.get_logger().info(f'  Description: {desc}')
        
        self.tf_pub = ros2_pat.utils.TransformPublisher(self)
        self.tf_pub.send_w_enu_to_ned_transform(rclpy.clock.Clock().now())

        self.vehicle_marker_pub = ros2_pat.utils.MarkerPublisher(self, '/pat/reference', mtype=10,
                                                                 mesh='package://ros2_pat/meshes/quad.dae',
                                                                 #mesh="package://ros2_pat/meshes/quadcopter.dae",
                                                                 #mesh="package://ros2_pat/meshes/mavic.dae",
                                                                 color=[0., 1., 0., 1.], frame_id="/w_ned")


        try:
            nwp = len(self.traj._g.waypoints)
            self.wp_pub = ros2_pat.utils.MarkerArrayPublisher(self, 'pat/waypoints', mtype=2,
                                                              meshes=['']*nwp,
                                                              colors=[[1., 0.2, 0.2, 0.5]]*nwp, scales=[(0.1, 0.1, 0.1)]*nwp)
        except AttributeError: pass
            
        self.trj_pub = ros2_pat.utils.TrajPublisher(self.traj, self, 'pat/ref_trajectory', color=[0.5, 0.75, 0.75, 0.25])
            
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
        #T_ned2b = ros2_pat.utils.T_of_t_q(X[fdm.sv_slice_pos], X[fdm.sv_slice_quat])
        # from fdm
        #sv_qi   = 6  # world to body rotation quaternion
        #sv_qx   = 7  #
        #sv_qy   = 8  #
        #sv_qz   = 9  #
        qalt = [X[fdm.sv_qx],X[fdm.sv_qy],X[fdm.sv_qz],X[fdm.sv_qi]]
        T_ned2b = ros2_pat.utils.T_of_t_q(X[fdm.sv_slice_pos], qalt)
        self.vehicle_marker_pub.publish(now, T_ned2b)
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
