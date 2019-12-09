#! /usr/bin/env python
import os, numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg
import dynamic_reconfigure.server
import pdb

import pat3.utils as p3_u
import pat3.algebra as p3_al
import pat3.frames as p3_fr
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn

# should be in ros_pat
import pat3.ros_utils as p3_rpu
import ros_pat.cfg.glider_trimConfig

class Agent(p3_rpu.PeriodicNode):
    def __init__(self, name='ros_pat_glider_trim'):
        p3_rpu.PeriodicNode.__init__(self, name)
        self.joint_state_pub = p3_rpu.JointStatePublisher(what=name)
        #self.ivel_pub = p3_rpu.PosePublisher(what=name)
        self.txt_pub = p3_rpu.TextMarkerPublisher('/ds_glider/trim_txt', 'ds_glider/base_link', scale=0.1, what=name)
        self.tf_pub = p3_rpu.TransformPublisher()
        self.tf_pub.send_w_enu_to_ned_transform(rospy.Time.now())

        ac_name = 'cularis'
        self.dm = p1_fw_dyn.DynamicModel(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))
        self.trim_args = {'h':0, 'va':10, 'throttle':0}
        self._trim()
        self.dyn_cfg_srv = dynamic_reconfigure.server.Server(ros_pat.cfg.glider_trimConfig, self.dyn_cfg_callback)

    def _trim(self):
        self.Xe, self.Ue = self.dm.trim(self.trim_args, report=True, debug=False)
        self.T_w2b = np.eye(4)
        phi, theta, psi = self.Xe[p3_fr.SixDOFAeroEuler.sv_slice_eul]
        p3_u.set_rot(self.T_w2b, p3_al.rmat_of_euler([phi, theta, psi]))
        self.T_b2w = np.linalg.inv(self.T_w2b)
        self.T_a2b = np.eye(4)
        va, self.alpha, beta = self.Xe[p3_fr.SixDOFAeroEuler.sv_slice_vaero]
        self.gamma = theta - self.alpha
        p3_u.set_rot(self.T_a2b, p3_fr.R_aero_to_body(self.alpha, beta))
        self.dail, self.dele, self.dflap = 0, self.Ue[self.dm.iv_de()], self.Ue[self.dm.iv_df()]
         
    def dyn_cfg_callback(self, config, level):
        self.trim_args['va']    = config['va']
        self.trim_args['flaps'] = np.deg2rad(config['flaps_deg'])
        self.trim_args['throttle']    = config['throttle']
        self._trim()
        return config 
     
    def periodic(self):
        now = rospy.Time.now()
        self.tf_pub.send_transform('w_ned', 'ds_glider/base_link', now, self.T_b2w)
        self.tf_pub.send_transform('ds_glider/base_link', 'ds_glider/aero', now, self.T_a2b)
        if 0:
            _v1 = 0.25*np.sin(1.5* now.to_sec())
            _v2 = 0.5*np.sin(0.5* now.to_sec())
            _v3 = 0.5*np.sin(2.* now.to_sec())
            joint_states = [_v1, -_v1, _v3, _v3, _v2, _v2]
        else:
            joint_states = [self.dail, -self.dail, self.dele, self.dele, self.dflap, self.dflap]
            
        self.joint_state_pub.publish(joint_states, now)
        #self.ivel_pub.publish(None)
        self.txt_pub.publish('va {:.1f}m/s\nalpha {:.1f} deg\ngamma {:.1f} deg'.format(self.trim_args['va'], np.rad2deg(self.alpha), np.rad2deg(self.gamma)))

if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    Agent().run(30)

