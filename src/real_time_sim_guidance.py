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

class Agent:
    def __init__(self):
      _fdm = fdm.FDM()

      _ctl_in = ctl.PosStep()
      #_traj = pmt.Circle([0, 0, -1], r=1.5, v=0.1)
      #_ctl_in = ctl.TrajRef(_traj, _fdm)

      _ctl = ctl.PosController(_fdm, _ctl_in)
      self.sim = pmu.Sim(_fdm, _ctl)
      self.tf_pub = pru.TransformPublisher()
      self.pose_pub = pru.QuadAndRefPublisher()
      
    def periodic(self):
        now = rospy.Time.now()
        self.sim.run(now.to_sec())
        self.tf_pub.publish(now, self.sim.fdm.T_w2b)
        self.pose_pub.publish([self.sim.fdm.T_w2b, self.sim.ctl.T_w2b_ref])
    
    def run(self):
        rate = rospy.Rate(20.)
        self.sim.reset(rospy.Time.now().to_sec())
        try:
            while not rospy.is_shutdown():
                self.periodic()
                rate.sleep()
        except rospy.exceptions.ROSInterruptException:
            pass  

def main():
    rospy.init_node('real_time_sim')
    Agent().run()
    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
