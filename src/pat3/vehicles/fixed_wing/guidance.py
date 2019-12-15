import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil

import pat3.utils as p3_u
import pat3.algebra as p3_alg
import pat3.frames as p3_fr
import pat3.trajectory_3D as p3_traj3d

class FMS:
    mod_auto1, mod_circle, mod_centering, mod_searching, mod_nb = range(5)
    def __init__(self, dm, trim_args, dt=0.01):
        self.mode = FMS.mod_circle
        self.Xe, self.Ue = dm.trim(trim_args, report=True)
        self.guidances = [ GuidanceAuto1(self.Xe, self.Ue, dm, dt),
                           GuidanceCircle(dm, trim_args, dt=0.01),
                           GuidanceThermal(dm, trim_args),
                           GuidanceSearching(dm, trim_args)]
        self.spv_size = 2
    
    def set_mode(self, m):
        if self.mode == FMS.mod_circle and m == FMS.mod_centering:
            #self.thermaling_traj.c = self.traj.c
            #self.thermaling_traj.r = self.traj.r
            self.guidances[FMS.mod_centering].set_circle_center(self.guidances[FMS.mod_circle].traj.c)
        self.mode = m   

    def guidance(self): return self.guidances[self.mode]
        
    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        return self.guidances[self.mode].get(t, X, Xee, Yc, debug, traj)

    def carrot(self): return self.guidance().carrot
    def phi_sp(self): return self.guidance().phi_sp
    def theta_sp(self): return self.guidance().theta_sp
    def get_va_sp(self): return self.guidance().v_sp
    def set_va_sp(self, _v): self.guidance().v_sp = _v


class GuidanceAuto1:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = Xe, Ue
        self.att_ctl = p3_pil.AttCtl(self.Xe, self.Ue, dm, dt)
        self.phi_sp = Xe[p3_fr.SixDOFAeroEuler.sv_phi]
        self.theta_sp = Xe[p3_fr.SixDOFAeroEuler.sv_theta]
        self.v_sp = self.Xe[p3_fr.SixDOFAeroEuler.sv_va]
        self.carrot= (0, 0, 0)          # unused here
        self.traj = p3_traj3d.RefTraj() # unused too

    def set_setpoints(self, phi_sp, theta_sp):
        #self.phi_sp, self.theta_sp = phi_sp, theta_sp
        self.phi_sp = phi_sp
        
    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        U = self.att_ctl.get(t, X, self.phi_sp, self.theta_sp)
        return U
     
        
class GuidancePurePursuit:
    v_mode_throttle, v_mode_vz, v_mode_alt = range(3)

    def __init__(self, dm, traj, trim_args={'h':0, 'va':12, 'gamma':0}, dt=0.01):
        self.Xe, self.Ue = dm.trim(trim_args, report=True)
        self.traj = traj
        self.carrot = np.zeros(3)
        self.phi_sp, self.theta_sp = self.Xe[p3_fr.SixDOFAeroEuler.sv_phi], self.Xe[p3_fr.SixDOFAeroEuler.sv_theta]
        self.att_ctl = p3_pil.AttCtl(self.Xe, self.Ue, dm, dt)
        self.T_w2b_ref = np.eye(4)
        self.sum_err_z, self.sum_err_v = 0, 0
        self.v_mode = self.v_mode_throttle#self.v_mode_alt
        self.v_sp = trim_args['va']
        self.throttle_sp = 0.
     
    def _compute_roll_setpoint(self, t, X, traj, max_phi=np.deg2rad(45)):
        my_pos = X[p3_fr.SixDOFAeroEuler.sv_slice_pos]
        self.carrot = traj.get_point_ahead(my_pos, 15.)
        self.b2c_ned = self.carrot-my_pos
        R_w2b = p3_alg.rmat_of_euler(X[p3_fr.SixDOFAeroEuler.sv_slice_eul])
        self.b2c_b = np.dot(R_w2b, self.b2c_ned)
        phi, theta = X[p3_fr.SixDOFAeroEuler.sv_phi], X[p3_fr.SixDOFAeroEuler.sv_theta]
        if 0:
            self.R = (np.linalg.norm(self.b2c_b)**2)/(2*self.b2c_b[1])
            # R.g.tan(phi) = v^2
            phi_sp = np.arctan(v**2/self.R/9.81)
        else:
            self.R = 0
            err_psi = X[p3_fr.SixDOFAeroEuler.sv_psi] - np.arctan2(self.b2c_ned[1], self.b2c_ned[0])
            err_psi = p3_u.norm_mpi_pi(err_psi)
            phi_sp = -0.75*err_psi
        self.phi_sp = np.clip(phi_sp, -max_phi, max_phi)
        return phi_sp
        
    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        if traj is None: traj = self.traj
        self._compute_roll_setpoint(t, X, traj)
        z, v = X[p3_fr.SixDOFAeroEuler.sv_z], X[p3_fr.SixDOFAeroEuler.sv_va]
        self.theta_sp = self.Xe[p3_fr.SixDOFAeroEuler.sv_theta] + 0.000025*self.sum_err_v
        U = self.att_ctl.get(t, X, self.phi_sp, self.theta_sp)
        _thr, _ail, _ele = 0, 1, 2
        # elevator compensation for banking
        v_sp, z_sp, theta_sp = self.v_sp, self.Xe[p3_fr.SixDOFAeroEuler.sv_z], self.Xe[p3_fr.SixDOFAeroEuler.sv_theta]
        self.sum_err_v += (v-v_sp)
        self.sum_err_z += (z-z_sp)
        d_ele = 0#self.sum_err_v*-0.00001 # -np.deg2rad(3.)*np.abs(self.phi_sp)/max_phi
        U[_ele] += d_ele
        if self.v_mode == self.v_mode_throttle:
            U[_thr] = self.throttle_sp
        else: # throttle integrator
            d_throttle = self.sum_err_z*0.000005
            U[_thr] += d_throttle
        if debug:
            fmt = 'z {:.1f} v {:.1f} m/s phi {:.1f}/{:.1f} deg theta {:.1f}/{:.1f} deg dthrottle {:.1f} % dele {:.1f} deg s_ev: {:.0f} s_ez: {:.0f}'
            print(fmt.format(z, v, np.rad2deg(phi), np.rad2deg(self.phi_sp), np.rad2deg(theta), np.rad2deg(self.theta_sp), d_throttle*100, np.rad2deg(d_ele), self.sum_err_v, self.sum_err_z))
        
        return U


class GuidanceCircle(GuidancePurePursuit):
    def __init__(self, dm, trim_args = {'h':0, 'va':12, 'gamma':0}, dt=0.01):
        traj = p3_traj3d.CircleRefTraj(c=[0., 0., 0.], r=15.)
        GuidancePurePursuit.__init__(self, dm, traj, trim_args, dt)

    def set_params(self, _c, _r):
        self.traj.c = np.asarray(_c)
        self.traj.r = _r


class GuidanceThermal(GuidancePurePursuit):

    def __init__(self, dm, traj, trim_args={'h':0, 'va':12, 'gamma':0}, dt=0.01):
        GuidancePurePursuit.__init__(self, dm, p3_traj3d.CircleRefTraj(c=[0., 0., 0.], r=15.), trim_args, dt)
        self.v_mode_throttle = GuidancePurePursuit.v_mode_throttle
        self.throttle_sp = 0.
        self.k_centering = 0.002
        # circle discretization
        self.nb_angles = 360
        self.avg_climb = np.zeros(self.nb_angles)
        self.avg_climb[:] = np.float('nan')

    def _idx_of_angle(self, angle): return int(p3_u.norm_0_2pi(angle)/2/np.pi*float(self.nb_angles))
    def _angle_of_idx(self, idx): return p3_u.norm_mpi_pi(idx/float(self.nb_angles)*2*np.pi)    

    def set_param(self, _v): self.k_centering = _v
    def set_radius(self, _v): self.traj.r = _v
    def set_circle_center(self, _c): self.traj.c = np.asarray(_c)
    
         
    def _thermal_centering(self, t, X, Xee):
        X_pos = X[p3_fr.SixDOFAeroEuler.sv_slice_pos]
        # climb rate
        ivel_ned = Xee[p3_fr.SixDOFEuclidianEuler.sv_slice_vel]
        alpha = self.traj.get_alpha(X_pos)
        #pdb.set_trace()
        idx_alpha = self._idx_of_angle(alpha)
        lp = 0.25
        #print(alpha, idx_alpha, ivel_ned[2])
        self.avg_climb[idx_alpha] = ivel_ned[2]#lp*self.avg_climb[idx_alpha] + (1-lp)*ivel_ned[2]
        #print(self.avg_climb[idx_alpha])
        if np.isnan(self.avg_climb).any():
            self.delta_center = np.zeros(3)
            #print('building gradiant')
        else:
            amin, amax = np.argmin(self.avg_climb), np.argmax(self.avg_climb)
            alpha_min, alpha_max = self._angle_of_idx(amin), self._angle_of_idx(amax)
            d_climb = self.avg_climb[amax] - self.avg_climb[amin]
            #fmt = 'alpha:{:.1f}deg min {:.0f}deg {:.2f}m/s max {:.0f}deg {:.2f}m/s diff: {:.1f}m/s'
            #args = np.rad2deg(alpha), np.rad2deg(alpha_min), self.avg_climb[amin], np.rad2deg(alpha_max), self.avg_climb[amax], d_climb
            #print(fmt.format(args))
            alpha_max_grad = alpha_max#(alpha_max+alpha_min + np.pi)/2
            self.delta_center = np.array([np.cos(alpha_max_grad), np.sin(alpha_max_grad), 0.])*d_climb*self.k_centering
            
    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        # offset circle toward vz gradient
        self._thermal_centering(t, X, Xee)
        self.traj.offset(self.delta_center)
        
        self.traj.c[2] = Xee[p3_fr.SixDOFEuclidianEuler.sv_z]
        return GuidancePurePursuit.get(self, t, X, Xee, Yc, debug)

    
class GuidanceSearching(GuidancePurePursuit):
    def __init__(self, dm, traj, trim_args = {'h':0, 'va':12, 'gamma':0}, dt=0.01):
        GuidancePurePursuit.__init__(self, dm, p3_traj3d.CircleRefTraj(c=[0., 0., 0.], r=50.), trim_args, dt)
