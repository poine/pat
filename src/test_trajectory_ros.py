#! /usr/bin/env python
import numpy as np
import rospy

import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl
import pat3.ros_utils as pru
import pat3.algebra as pal

class Agent:

    def __init__(self, traj, time_factor=1.):
        self.time_factor = time_factor
        self.traj = traj
        self.tf_pub = pru.TransformPublisher()
        self.pose_pub = pru.PoseArrayPublisher()
        self.poles_pub = pru.TrackPublisher()
        self.traj_pub = pru.TrajectoryPublisher(self.traj)
        self.fdm = fdm.FDM()
        self.df = ctl.DiffFlatness()
        
    def periodic(self):
        now = rospy.Time.now()
        #print('periodic at {}'.format(now.to_sec()))
        t_sim = self.t0 + (now.to_sec() - self.t0)*self.time_factor
        Yc = self.traj.get(t_sim)
        Xc, Uc = self.df.state_and_cmd_of_flat_output(Yc, self.fdm.P)
        
        T_w2b = np.eye(4)
        T_w2b[:3,:3] = pal.rmat_of_quat(Xc[fdm.sv_slice_quat]).T # that is freaking weird....
        T_w2b[:3,3] = Yc[:3,0]
        self.pose_pub.publish([T_w2b])
        self.tf_pub.publish(now, T_w2b)
        self.poles_pub.publish()
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
    #traj = pmt.Circle([0, 0, 0], 1, v=4.)
    #traj = pmt.Line([0, 0, 0], [10, 0, 0], 0.1)
    #traj = pmt.Oval(l=1, r=1, v=5)
    #traj = pmt.DoubleOval(l=2, r=1, v=5)
    if 0:
        r, v = 1., 4.; om = v/r
        traj = pmt.Circle(c=[0, 0, 1], r=1., v=4., alpha0=0., dalpha=2*np.pi, zt=pmt.SinOne(a=1.5, om=om))
    #traj = pmt.Circle(c=[0, 0, -1], r=1., v=4., alpha0=0., dalpha=2*np.pi)
    #traj = pmt.Circle(c=[0, 0, -1], r=1., v=4., alpha0=0., dalpha=2*np.pi, psit=pmt.CstOne(0))
    traj =  pmt.FigureOfEight(v=1.)
    Agent(traj, time_factor=0.5).run()
    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
