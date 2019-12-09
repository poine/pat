#! /usr/bin/env python
import numpy as np
import rospy, tf2_ros, tf, geometry_msgs.msg, sensor_msgs.msg

import pdb

import pat3.utils as p3_u
import pat3.algebra as p3_al
import pat3.frames as p3_fr
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
# should be in ros_pat
import pat3.ros_utils as p3_rpu

import pat3.test.fixed_wing.test_01_dynamics as test_dyn
import pat3.test.fixed_wing.test_02_att_ctl  as test_pil


class Agent(p3_rpu.PeriodicNode):
    mode_cst, mode_sim = range(2)

    def __init__(self, mode=mode_sim):
        p3_rpu.PeriodicNode.__init__(self, 'ros_pat_play_simulation')
        self.mode = mode

        pos_ned = [0, 0, -1]
        self.T_w2b = np.eye(4)
        self.T_b2w = np.linalg.inv(self.T_w2b)
        self.T_a2b = np.eye(4)
        self.T_b2a = self.T_a2b#np.linalg.inv(self.T_a2b) # FIXME
        save=False
        #savefile_name = '/tmp/pat_traj.npz'
        savefile_name = '/tmp/pat_glider_circle.npz'
        if mode == self.mode_sim:
            if save:
                #self.time, self.X = test_dyn.sim_step_thr()
                #self.time, self.X = test_dyn.sim_step_ail()
                #self.time, self.X = test_dyn.sim_step_ele()
                self.time, self.X = test_dyn.sim_step_rud()
                #self.time, self.X, self.U = test_pil.test_step_theta(dm)
                np.savez(savefile_name, time=self.time, X=self.X)
                print('saved {}'.format(savefile_name))
            else:
                print('loading {}'.format(savefile_name))
                _data =  np.load(savefile_name)
                self.time, self.X =[_data[k] for k in ['time', 'X']] 
            self.sim_dt = self.time[1] - self.time[0]
            self.sim_dur = self.time[-1] - self.time[0]
            
        self.tf_pub = p3_rpu.TransformPublisher()
        self.marker_pub = p3_rpu.PoseArrayPublisher(dae='ds_glider_full.dae')
        self.odom_pub = p3_rpu.OdomPublisher()
        
        
    def periodic(self):
        now = rospy.Time.now()
        self.tf_pub.send_w_enu_to_ned_transform(now)
        if self.mode == self.mode_sim:
            t_sim = np.fmod(now.to_sec(), self.sim_dur)
            idx_t = int(t_sim / self.sim_dt)
            X_eul = self.X[idx_t, p3_fr.SixDOFAeroEuler.sv_slice_eul]
            #X_eul[0] = np.sin(t_sim)
            #print t_sim, idx_t, self.X[idx_t, p3_fr.SixDOFAeroEuler.sv_slice_eul]
            p3_u.set_rot(self.T_w2b, p3_al.rmat_of_euler(X_eul))
            #_set_trans(self.T_w2b, self.X[idx_t, p3_fr.SixDOFAeroEuler.sv_slice_pos])
            self.T_b2w = np.linalg.inv(self.T_w2b)
            p3_u.set_trans(self.T_b2w, self.X[idx_t, p3_fr.SixDOFAeroEuler.sv_slice_pos])
            va, alpha, beta = self.X[idx_t, p3_fr.SixDOFAeroEuler.sv_slice_vaero]
            p3_u.set_rot(self.T_a2b, p3_fr.R_aero_to_body(alpha, beta))
            self.T_b2a = self.T_a2b#np.linalg.inv(self.T_a2b) # FIXME
            
                
        self.tf_pub.send_w_ned_to_b_transform(now, self.T_b2w)
        self.tf_pub.send_b_to_a_transform(now, self.T_b2a)
        self.tf_pub.send_transform('w_ned', 'ds_glider/base_link', now, self.T_b2w)
        self.marker_pub.publish([self.T_b2w])
        self.odom_pub.publish(self.T_b2w, now)
        
def main():
    Agent().run(25)


if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
