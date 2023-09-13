import rclpy
import rclpy.node
import numpy as np

import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as trjf
import ros2_pat.utils

import pdb

class WpNavPublisher(rclpy.node.Node):
    def __init__(self, wp_traj):
        super().__init__('WpNavPublisher')
        self.wp_traj = wp_traj
        print(wp_traj._g.waypoints)
        nwp = len(wp_traj._g.waypoints)
        self.wp_pub = ros2_pat.utils.MarkerArrayPublisher(self, 'pat/waypoints', mtype=2, meshes=['']*nwp, colors=[[1., 0.2, 0.2, 0.5]]*nwp, scales=[(0.1, 0.1, 0.1)]*nwp)

        self.trj_pub = ros2_pat.utils.TrajPublisher(wp_traj, self, 'pat/ref_trajectory')
        

        timer_period = 1./2.  # seconds
        self.timer = self.create_timer(timer_period, self.timer_callback)
        
        
    def timer_callback(self):
        nwp = len(self.wp_traj._g.waypoints)
        T_ned2b = np.array([np.eye(4)]*nwp)
        for i in range(nwp):
            T_ned2b[i][:3,3] = self.wp_traj._g.waypoints[i]
            #print(self.wp_traj._g.waypoints[i], T_ned2b[i])
        self.wp_pub.publish(T_ned2b)
        self.trj_pub.publish()
        
        
    
def main(args=None):
    rclpy.init(args=args)
    wp_publisher = WpNavPublisher(trjf.Traj43())
    rclpy.spin(wp_publisher)
    wp_publisher.destroy_node()
    rclpy.shutdown()


if __name__ == '__main__':
    main()
