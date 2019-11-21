#! /usr/bin/env python
import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil

import pat3.atmosphere as p3_atm
import pat3.frames as p3_fr
import pat3.trajectory_3D as p3_traj3d
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.algebra as p3_alg

import control.matlab

import test_02_att_ctl


#def norm_mpi_pi(v): return ( v + np.pi) % (2 * np.pi ) - np.pi

class Guidance:
    v_mode_throttle, v_mode_vz, v_mode_alt = range(3)

    def __init__(self, dm, traj, trim_args = {'h':0, 'va':12, 'gamma':0}, dt=0.01):
        self.Xe, self.Ue = dm.trim(trim_args, debug=True)
        self.traj = traj
        self.att_ctl = p3_pil.AttCtl(self.Xe, self.Ue, dm, dt)
        self.spv_size = 2
        self.T_w2b_ref = np.eye(4)
        self.sum_err_z, self.sum_err_v = 0, 0
        self.v_mode = self.v_mode_throttle#self.v_mode_alt
        self.v_sp = 12
 

    def _idx_of_angle(self, angle): return int(p3_u.norm_0_2pi(angle)/2/np.pi*float(self.nb_angles))
    def _angle_of_idx(self, idx): return p3_u.norm_mpi_pi(idx/float(self.nb_angles)*2*np.pi)
    
    def set_circle(self, xc, yc, zc, r):
        self.traj.c = np.array([xc, yc, zc])
        self.traj.r = r

    def _compute_roll_setpoint(self, t, X, traj, max_phi=np.deg2rad(45)):
        my_pos = X[p1_fw_dyn.sv_slice_pos]
        self.carrot = traj.get_point_ahead(my_pos, 15.)
        self.b2c_ned = self.carrot-my_pos
        R_w2b = p3_alg.rmat_of_euler(X[p1_fw_dyn.sv_slice_eul])
        self.b2c_b = np.dot(R_w2b, self.b2c_ned)
        phi, theta = X[p1_fw_dyn.sv_phi], X[p1_fw_dyn.sv_theta]
        if 0:
            self.R = (np.linalg.norm(self.b2c_b)**2)/(2*self.b2c_b[1])
            # R.g.tan(phi) = v^2
            phi_sp = np.arctan(v**2/self.R/9.81)
        else:
            self.R = 0
            err_psi = X[p1_fw_dyn.sv_psi] - np.arctan2(self.b2c_ned[1], self.b2c_ned[0])
            err_psi = p3_u.norm_mpi_pi(err_psi)
            phi_sp = -0.75*err_psi
        self.phi_sp = np.clip(phi_sp, -max_phi, max_phi)
        return phi_sp

            
        
        
    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        if traj is None: traj = self.traj
        phi_sp = self._compute_roll_setpoint(t, X, traj)
        z, v = X[p1_fw_dyn.sv_z], X[p1_fw_dyn.sv_v]
        self.theta_sp = self.Xe[p1_fw_dyn.sv_theta] + 0.000025*self.sum_err_v
        U = self.att_ctl.get(t, X, self.phi_sp, self.theta_sp)
        _thr, _ail, _ele = 0, 1, 2
        # elevator compensation for banking
        v_sp, z_sp, theta_sp = self.v_sp, self.Xe[p1_fw_dyn.sv_z], self.Xe[p1_fw_dyn.sv_theta]
        self.sum_err_v += (v-v_sp)
        self.sum_err_z += (z-z_sp)
        d_ele = 0#self.sum_err_v*-0.00001 # -np.deg2rad(3.)*np.abs(phi_sp)/max_phi
        U[_ele] += d_ele
        if self.v_mode == self.v_mode_throttle:
            U[_thr] = 0.0
        else: # throttle integrator
            d_throttle = self.sum_err_z*0.000005
            U[_thr] += d_throttle
        if debug:
            fmt = 'z {:.1f} v {:.1f} m/s phi {:.1f}/{:.1f} deg theta {:.1f}/{:.1f} deg dthrottle {:.1f} % dele {:.1f} deg s_ev: {:.0f} s_ez: {:.0f}'
            print(fmt.format(z, v, np.rad2deg(phi), np.rad2deg(self.phi_sp), np.rad2deg(theta), np.rad2deg(self.theta_sp), d_throttle*100, np.rad2deg(d_ele), self.sum_err_v, self.sum_err_z))
        
        return U


class GuidanceThermal(Guidance):
    fms_mode_circle, fms_mode_thermaling = range(2)
    def __init__(self, dm, traj, trim_args = {'h':0, 'va':12, 'gamma':0}, dt=0.01):
        Guidance.__init__(self, dm, traj, trim_args, dt)
        self.v_mode_throttle = Guidance.v_mode_throttle
        self.throttle_sp = 0.
        self.thermaling_traj = p3_traj3d.CircleRefTraj(c=[0., 0., 0.], r=15.)
        self.fms_mode = GuidanceThermal.fms_mode_circle
        self.k_centering = 0.002
        # circle discretization
        self.nb_angles = 360
        self.avg_climb = np.zeros(self.nb_angles)
        self.avg_climb[:] = np.float('nan')
        
    def set_mode(self, m):
        if self.fms_mode == GuidanceThermal.fms_mode_circle and m == GuidanceThermal.fms_mode_thermaling:
            self.thermaling_traj.c = self.traj.c
            self.thermaling_traj.r = self.traj.r
        self.fms_mode = m
    def set_param(self, _v): self.k_centering = _v
    def set_radius(self, _v): self.thermaling_traj.r = _v
        
    def _thermal_centering(self, t, X, Xee):
        X_pos = X[p1_fw_dyn.sv_slice_pos]
        # climb rate
        ivel_ned = Xee[p3_fr.SixDOFEuclidianEuler.sv_slice_v]
        alpha = self.traj.get_alpha(X_pos)
        #pdb.set_trace()
        idx_alpha = self._idx_of_angle(alpha)
        lp = 0.25
        #print(alpha, idx_alpha, ivel_ned[2])
        self.avg_climb[idx_alpha] = ivel_ned[2]#lp*self.avg_climb[idx_alpha] + (1-lp)*ivel_ned[2]
        #print(self.avg_climb[idx_alpha])
        amin, amax = np.argmin(self.avg_climb), np.argmax(self.avg_climb)
        alpha_min, alpha_max = self._angle_of_idx(amin), self._angle_of_idx(amax)
        d_climb = self.avg_climb[amax] - self.avg_climb[amin]
        #fmt = 'alpha:{:.1f}deg min {:.0f}deg {:.2f}m/s max {:.0f}deg {:.2f}m/s diff: {:.1f}m/s'
        #args = np.rad2deg(alpha), np.rad2deg(alpha_min), self.avg_climb[amin], np.rad2deg(alpha_max), self.avg_climb[amax], d_climb
        #print(fmt.format(args))
        alpha_max_grad = alpha_max#(alpha_max+alpha_min + np.pi)/2
        self.delta_center = np.array([np.cos(alpha_max_grad), np.sin(alpha_max_grad), 0.])*d_climb*self.k_centering
        #if t > 20:
            #print('{} {}'.format(delta_center, self.traj.c))#np.rad2deg(alpha_max_grad)))

        
    def get(self, t, X, Xee=None, Yc=None, debug=False):
        try:
            self._thermal_centering(t, X, Xee)
            if self.fms_mode == GuidanceThermal.fms_mode_thermaling:
                self.thermaling_traj.offset(self.delta_center)
                self.thermaling_traj.c[2] = Xee[p3_fr.SixDOFEuclidianEuler.sv_z]
                return Guidance.get(self, t, X, Xee, Yc, debug, self.thermaling_traj)
        except ValueError:
            print('not enought data to compute gradiant')
        # default to circle
        return Guidance.get(self, t, X, Xee, Yc, debug)
    
    
def run_simulation(dm, ref_traj, tf=60.5, dt=0.01, trim_args = {'h':0, 'va':12, 'gamma':0}, plot=False, atm=None):
    time = np.arange(0, tf, dt)
    X = np.zeros((len(time), dm.sv_size))
    U = np.zeros((len(time),  dm.input_nb()))
    carrots, att_sp = np.zeros((len(time),  3)), np.zeros((len(time), 2))
    ctl = Guidance(dm, ref_traj, trim_args, dt)
    X[0] = dm.reset(ctl.Xe, t0=time[0], X_act0=None)#ctl.Ue)
    for i in range(1, len(time)):
        Xee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(X[i-1])
        U[i-1] = ctl.get(time[i-1], X[i-1], Xee)
        carrots[i-1] = ctl.carrot; att_sp[i-1] = ctl.phi_sp, ctl.theta_sp 
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1], atm)
    U[-1] = ctl.get(time[-1], X[-1]); carrots[-1] = ctl.carrot; att_sp[-1] = ctl.phi_sp, ctl.theta_sp 
    if plot:
        dm.plot_trajectory(time, X, U)
        plt.subplot(5,3,1); plt.plot(time, carrots[:,0])
        plt.subplot(5,3,2); plt.plot(time, carrots[:,1])
        plt.subplot(5,3,3); plt.plot(time, carrots[:,2])
        plt.subplot(5,3,7); plt.plot(time, np.rad2deg(att_sp[:,0]), label='setpoint')
        plt.subplot(5,3,8); plt.plot(time, np.rad2deg(att_sp[:,1]), label='setpoint')
        plt.show()  
    return time, X, U


    
    
def main(param_filename, trim_args = {'h':0, 'va':11, 'gamma':0}):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    ref_traj = p3_traj3d.CircleRefTraj(c=[0, 0, 0], r=20)
    #time, X, U = test_02_att_ctl.run_simulation(dm, plot=False)
    time, X, U = run_simulation(dm, ref_traj, plot=True)
    
    #atm =  p3_atm.Atmosphere()
    p3_pu.plot_3D_traj(ref_traj, X)
    #p3_pu.plot_3D_wind(atm)
    plt.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    dirname, filename = os.path.split(os.path.abspath(__file__))
    param_filename = os.path.abspath(os.path.join(dirname, '../../../../data/vehicles/cularis.xml'))
    main(param_filename)
