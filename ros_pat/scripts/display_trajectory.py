#! /usr/bin/env python
import numpy as np
import rospy

import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as pmtf
#import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl
import pat3.algebra as pal

import pat3.ros_utils as pru

'''
This is a ROS agent that sends messages for displaying a trajectory in rviz
'''

class Agent:

    def __init__(self, traj, time_factor=1.):
        self.time_factor = time_factor
        self.traj = traj
        self.tf_pub = pru.TransformPublisher()
        #self.pose_pub = pru.PoseArrayPublisher(dae='quad.dae')
        self.pose_pub = pru.PoseArrayPublisher(dae='glider.dae')
        self.traj_pub = pru.TrajectoryPublisher(self.traj)
        self.fdm = fdm.FDM()
        self.df = ctl.DiffFlatness()
        
    def periodic(self):
        now = rospy.Time.now()
        #print('periodic at {}'.format(now.to_sec()))
        t_sim = self.t0 + (now.to_sec() - self.t0)*self.time_factor
        Yc = self.traj.get(t_sim)
        Xc, Uc, Xdc = self.df.state_and_cmd_of_flat_output(Yc, self.fdm.P)
        
        T_w2b = np.eye(4)
        T_w2b[:3,:3] = pal.rmat_of_quat(Xc[fdm.sv_slice_quat]).T # that is freaking weird....
        T_w2b[:3,3] = Yc[:3,0]
        self.pose_pub.publish([T_w2b])
        self.tf_pub.publish(now, T_w2b)
        self.traj_pub.publish()                
    
    def run(self):
        rate = rospy.Rate(40.)
        self.t0 = rospy.Time.now().to_sec()
        self.traj.reset(self.t0)
        try:
            while not rospy.is_shutdown():
                self.periodic()
                rate.sleep()
        except rospy.exceptions.ROSInterruptException:
            pass

def main():
    rospy.init_node('test_trajectory')
    traj_name = rospy.get_param('~traj_name', 'circle4')
    if rospy.get_param('~list_available', False):
        rospy.loginfo('  available trajectories\n{}'.format(pmtf.list()))
    try:
        rospy.loginfo('  Loading trajectory: {}'.format(traj_name))
        traj, desc = pmtf.get(traj_name)
        rospy.loginfo('  Description: {}'.format(desc))
        Agent(traj, time_factor=0.5).run()
    except KeyError:
        rospy.loginfo('  Unkwown trajectory: {}'.format(traj_name))
        rospy.loginfo( '  available trajectories\n{}'.format(pmtf.list()))
        
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
