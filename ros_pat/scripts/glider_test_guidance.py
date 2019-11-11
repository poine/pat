#! /usr/bin/env python
import numpy as np
import rospy

import pat.vehicles.fixed_wing.dynamic_model_python_basic as p1_fw_dyn
import pat3.test.fixed_wing.test_03_guidance as p3_fw_guid
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
        traj = p3_fw_guid.CircleRefTraj2(c=[0, 0, 10], r=15)
        self.time_factor = time_factor
        self.traj = traj
        self.tf_pub = pru.TransformPublisher()
        self.pose_pub = pru.PoseArrayPublisher(dae='glider.dae')
        self.traj_pub = pru.TrajectoryPublisher(self.traj)

        self.carrot_pub  = pru.MarkerPublisher('guidance/goal', 'w_ned')
        self.carrot_pub2 = pru.MarkerPublisher('guidance/goal2', 'b_frd', argb=(1., 1., 1., 0.))
 
        self.fdm = fdm.FDM()
        param_filename='/home/poine/work/pat/data/vehicles/cularis.xml'
        self.dm = p1_fw_dyn.DynamicModel(param_filename)
        self.guidance = p3_fw_guid.Guidance(self.dm, traj, {'h':0, 'va':15, 'gamma':0})
        
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

        X1 = np.zeros(p1_fw_dyn.sv_size)
        X1[p1_fw_dyn.sv_slice_pos] = Xc[fdm.sv_slice_pos]
        X1[p1_fw_dyn.sv_slice_eul] = pal.euler_of_quat(Xc[fdm.sv_slice_quat])
        U = self.guidance.get(t_sim, X1)
        self.carrot_pub.publish(self.guidance.carrot)
        #b2c_bfrd = np.dot(T_w2b[:3, :3].T, self.guidance.b2c_ned)
        #print self.guidance.carrot
        #b2c_bfrd = np.dot(np.linalg.inv(T_w2b), [self.guidance.carrot[0], self.guidance.carrot[1], self.guidance.carrot[2], 1])
        b2c_bfrd = self.guidance.b2c_b
        self.carrot_pub2.publish(b2c_bfrd)
        
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
        Agent(traj=None, time_factor=0.5).run()
    except KeyError:
        rospy.loginfo('  Unkwown trajectory: {}'.format(traj_name))
        rospy.loginfo( '  available trajectories\n{}'.format(pmtf.list()))
        
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
