#! /usr/bin/env python
import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg

import pdb

import pat3.algebra as pal
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl


class Sim:

    def __init__(self):
        self.dt_internal = 0.01
        self.fdm = fdm.FDM()
        self.ctl = ctl.FooController(self.fdm)
        self.Yc = np.zeros(5); self.Yc[1] = 1

    def reset(self, t):
        self.t = t
        self.Xe, self.Ue = self.fdm.trim()
        self.X0 = np.array(self.Xe)
        #self.X0[fdm.sv_zd] += 0.05
        #self.X0[fdm.sv_r] += 0.05
        #phi0, theta0, psi0 = 0, -np.deg2rad(0.1), 0.
        #self.X0[fdm.sv_qx:fdm.sv_qw+1] = pal.quat_of_euler([phi0, theta0, psi0])
        return self.fdm.reset(self.X0, self.t)
        
    def run(self, t):
        U = self.ctl.get(self.fdm.X, self.Yc)
        self.fdm.run(t, U)
        return U, self.fdm.X
        

class Agent:

    def __init__(self):
        self.sim = Sim()
        self.tfBcaster = tf2_ros.TransformBroadcaster()
        self.tfBuffer  = tf2_ros.Buffer()
        self.tfLstener = tf2_ros.TransformListener(self.tfBuffer)
        rospy.Subscriber('/joy', sensor_msgs.msg.Joy, self.on_joy_msg, queue_size=1)
        #self.input = ctl.StepZInput()
        self.input = ctl.StepEulerInput(pal.e_phi, _a=np.deg2rad(1.), p=4, dt=1)
        #self.input = ctl.StepEulerInput(pal.e_theta)
        #self.input = ctl.StepEulerInput(pal.e_psi, _a=np.deg2rad(30.), p=2, dt=1)
        #self.input = ctl.CstInput(0, [np.deg2rad(0.1), 0, 0])
        
    def send_w_enu_to_ned_transform(self, t):
        R_enu2ned = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
        T_enu2ned = np.eye(4); T_enu2ned[:3,:3] = R_enu2ned
        self.send_transform("w_enu", "w_ned", t, T_enu2ned)

    def send_w_ned_to_b_transform(self, t):
        self.send_transform("w_ned", "b_frd", t, self.sim.fdm.T_w2b)

    def send_transform(self, f1, f2, t, T_f1tof2):
        tf_msg = geometry_msgs.msg.TransformStamped()
        tf_msg.header.frame_id = f1
        tf_msg.child_frame_id = f2
        tf_msg.header.stamp = t
        _r = tf_msg.transform.rotation
        _r.x, _r.y, _r.z, _r.w = tf.transformations.quaternion_from_matrix(T_f1tof2)
        _t = tf_msg.transform.translation
        _t.x, _t.y, _t.z = T_f1tof2[:3,3]
        self.tfBcaster.sendTransform(tf_msg)

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
        self.send_w_enu_to_ned_transform(now)
        self.send_w_ned_to_b_transform(now)
    
    def run(self):
        rate = rospy.Rate(100.)
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
