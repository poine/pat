#! /usr/bin/env python
import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg

import pdb

import pat3.algebra as pal
import pat3.utils as pmu
import pat3.ros_utils as pru
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl
import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt
import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as pmtf


class Agent:
    def __init__(self, _traj=None, time_factor=1.):
        self.time_factor = time_factor
        _fdm = fdm.MR_FDM()
        #_fdm = fdm.UFOFDM()
        #_fdm = fdm.SolidFDM()
        self.traj_pub = pru.TrajectoryPublisher(_traj)
        _ctl_in = ctl.TrajRef(_traj, _fdm)
        _ctl = ctl.PosController(_fdm, _ctl_in)
        self.sim = pmu.Sim(_fdm, _ctl)
        self.tf_pub = pru.TransformPublisher()
        self.pose_pub = pru.QuadAndRefPublisher()
      
    def periodic(self):
        now = rospy.Time.now()
        t_sim = self.t0 + (now.to_sec() - self.t0)*self.time_factor
        self.sim.run(t_sim)
        self.tf_pub.publish(now, self.sim.fdm.T_w2b)
        self.pose_pub.publish([self.sim.fdm.T_w2b, self.sim.ctl.T_w2b_ref])
        if self.traj_pub is not None: self.traj_pub.publish() 
        
    def run(self):
        rate = rospy.Rate(20.)
        self.t0 = rospy.Time.now().to_sec()
        self.sim.reset(self.t0, self.sim.ctl.setpoint.get(self.t0)[2])
        try:
            while not rospy.is_shutdown():
                self.periodic()
                rate.sleep()
        except rospy.exceptions.ROSInterruptException:
            pass  

def main():
    rospy.init_node('real_time_sim')
    traj_name = rospy.get_param('~traj_name', 'circle_with_intro_slow')
    rospy.loginfo('  Loading trajectory: {}'.format(traj_name))
    traj, desc = pmtf.get(traj_name)
    time_factor = rospy.get_param('~time_factor', 1.)
    Agent(traj, time_factor).run()
    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
