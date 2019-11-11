#! /usr/bin/env python
import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg

import pdb

import pat3.algebra as pal
import pat3.utils as pmu
import pat3.ros_utils as pru
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl


class Agent:

    def __init__(self):
        _fdm = fdm.FDM()
        _ctl = ctl.ZAttController(_fdm)
        self.sim = pmu.Sim(_fdm, _ctl)
        self.tf_pub = pru.TransformPublisher()
        self.pose_pub = pru.PoseArrayPublisher()
        rospy.Subscriber('/joy', sensor_msgs.msg.Joy, self.on_joy_msg, queue_size=1)
        #self.input = ctl.StepZInput(0.1)
        #self.input = ctl.StepEulerInput(pal.e_phi, _a=np.deg2rad(1.), p=4, dt=1)
        self.input = ctl.StepEulerInput(pal.e_theta, _a=np.deg2rad(1.), p=4, dt=1)
        #self.input = ctl.StepEulerInput(pal.e_psi, _a=np.deg2rad(30.), p=2, dt=1)
        #self.input = ctl.CstInput(0, [np.deg2rad(0.1), 0, 0])
        #self.input = ctl.RandomInput()
        
 

    def on_joy_msg(self, msg):
        thr_val = msg.axes[2]
        pitch_val = msg.axes[1]
        roll_val = msg.axes[0]
        yaw_val = msg.axes[3]
        mod_val = msg.axes[7]
        #pdb.set_trace()
        #print(msg)
        self.sim.Yc[0] = thr_val
        
    def periodic(self):
        now = rospy.Time.now()
        self.sim.Yc = self.input.get(now.to_sec())
        self.sim.run(now.to_sec())
        self.tf_pub.publish(now, self.sim.fdm.T_w2b)
        self.pose_pub.publish([self.sim.fdm.T_w2b])
    
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
