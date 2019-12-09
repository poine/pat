#! /usr/bin/env python
import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg

import pdb

import pat3.utils as p3_u
import pat3.algebra as p3_al
import pat3.frames as p3_fr
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn

import pat3.ros_utils as p3_rpu
import pat3.test.fixed_wing.test_01_dynamics as test_dyn


class Agent(p3_rpu.PeriodicNode):

    def __init__(self):
        p3_rpu.PeriodicNode.__init__(self, 'ros_pat_publish_transform')
        pos_ned = [0, 0, -1]
        self.T_w2b = np.eye(4)
        dm_ae = p1_fw_dyn.DynamicModel()
        Xe, Ue = dm_ae.trim({'h':0, 'va':8, 'gamma':0}, report=True, debug=True)
        va, alpha, beta = Xe[p3_fr.SixDOFAeroEuler.sv_slice_vaero]
        #va, alpha, beta = 10, np.deg2rad(3), np.deg2rad(0)
        phi, theta, psi = Xe[p3_fr.SixDOFAeroEuler.sv_slice_eul]
        #phi, theta, psi = np.deg2rad(0), np.deg2rad(20), np.deg2rad(0) 
        p3_u.set_rot(self.T_w2b, p3_al.rmat_of_euler([phi, theta, psi]))
        self.T_b2w = np.linalg.inv(self.T_w2b)
        p3_u.set_trans(self.T_b2w, pos_ned)
        self.T_a2b = np.eye(4)
        p3_u.set_rot(self.T_a2b, p3_fr.R_aero_to_body(alpha, beta))
        self.T_b2a = self.T_a2b#np.linalg.inv(self.T_a2b) # FIXME

        self.tf_pub = p3_rpu.TransformPublisher()
        self.marker_pub = p3_rpu.PoseArrayPublisher(dae='glider.dae')
        
        
    def periodic(self):
        now = rospy.Time.now()
        self.tf_pub.send_w_enu_to_ned_transform(now)
        self.tf_pub.send_w_ned_to_b_transform(now, self.T_b2w)
        self.tf_pub.send_b_to_a_transform(now, self.T_b2a)
        self.tf_pub.send_transform('w_ned', 'fox/base_link', now, self.T_b2w)
        
        self.marker_pub.publish([self.T_b2w])
        
        
def main():
    Agent().run(10)


if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
