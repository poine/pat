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

      if 0:
          _ctl_in = ctl.PosStep()
          self.traj_pub = None
      else:
          #_traj = pmt.Circle([0, 0, -0.5], r=1.5, v=2.)
          _traj =  pmt.FigureOfHeight(v=1.)
          #_traj = pmt.Oval(l=1, r=1, v=2.)
          self.traj_pub = pru.TrajectoryPublisher(_traj)
          _ctl_in = ctl.TrajRef(_traj, _fdm)

      _ctl = ctl.PosController(_fdm, _ctl_in)
      self.sim = pmu.Sim(_fdm, _ctl)
      self.tf_pub = pru.TransformPublisher()
      self.pose_pub = pru.QuadAndRefPublisher()
      
    def periodic(self):
        now = rospy.Time.now()
        self.sim.run(now.to_sec())
        self.tf_pub.publish(now, self.sim.fdm.T_w2b)
        self.pose_pub.publish([self.sim.fdm.T_w2b, self.sim.ctl.T_w2b_ref])
        if self.traj_pub is not None: self.traj_pub.publish() 
        
    def run(self):
        rate = rospy.Rate(20.)
        t0 = rospy.Time.now().to_sec()
        self.sim.reset(t0, self.sim.ctl.setpoint.get(t0)[2])
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
