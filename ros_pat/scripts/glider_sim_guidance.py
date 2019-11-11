#! /usr/bin/env python
import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg

import pdb

import pat3.algebra as pal
import pat3.utils as pmu
import pat3.ros_utils as pru
import pat.vehicles.fixed_wing.dynamic_model_python_basic as p1_fw_dyn
import pat3.test.fixed_wing.test_03_guidance as p3_fw_guid
#import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
#import pat3.vehicles.rotorcraft.multirotor_control as ctl
import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt
import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as pmtf


class Agent:
    def __init__(self, _traj=None, time_factor=1.):
        traj = p3_fw_guid.CircleRefTraj2(c=[0, 0, 10], r=30)
        self.time_factor = time_factor
        #_fdm = fdm.MR_FDM()
        #_fdm = fdm.UFOFDM()
        #_fdm = fdm.SolidFDM()
        param_filename='/home/poine/work/pat/data/vehicles/cularis.xml'
        _fdm = p1_fw_dyn.DynamicModel(param_filename)
        self.traj_pub = pru.TrajectoryPublisher(_traj, ms=0.2)
        _ctl = p3_fw_guid.Guidance(_fdm, traj, {'h':0, 'va':12, 'gamma':0})
        self.sim = pmu.Sim(_fdm, _ctl)
        self.tf_pub = pru.TransformPublisher()
        self.carrot_pub = pru.MarkerPublisher('guidance/goal', 'w_ned', scale=(0.5, 0.5, 0.5))
        #self.pose_pub = pru.QuadAndRefPublisher()
        self.pose_pub = pru.PoseArrayPublisher(dae='glider.dae')

    def periodic(self):
        now = rospy.Time.now()
        t_sim = (now.to_sec() - self.t0)*self.time_factor
        self.sim.run(t_sim)
        self.tf_pub.publish(now, self.sim.fdm.T_w2b)
        #self.pose_pub.publish([self.sim.fdm.T_w2b, self.sim.ctl.T_w2b_ref])
        self.pose_pub.publish([self.sim.fdm.T_w2b])
        self.carrot_pub.publish(self.sim.ctl.carrot)
        if self.traj_pub is not None: self.traj_pub.publish() 
        
    def run(self):
        rate = rospy.Rate(20.)
        self.t0 = rospy.Time.now().to_sec()
        self.sim.reset(0, self.sim.ctl.Xe)
        try:
            while not rospy.is_shutdown():
                self.periodic()
                rate.sleep()
        except rospy.exceptions.ROSInterruptException:
            pass  

def main():
    rospy.init_node('real_time_sim')
    traj_name = rospy.get_param('~traj_name', 'c42')
    if rospy.get_param('~list_available', False):
        rospy.loginfo('  available trajectories\n{}'.format(pmtf.list()))
    rospy.loginfo('  Loading trajectory: {}'.format(traj_name))
    traj, desc = pmtf.get(traj_name)
    time_factor = rospy.get_param('~time_factor', 1.)
    rospy.loginfo('  Running simulation with time factor: {}'.format(time_factor))
    Agent(traj, time_factor).run()
    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
