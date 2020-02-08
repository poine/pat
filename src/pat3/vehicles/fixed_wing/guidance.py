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
import pat3.plot_utils as p3_pu



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

    def enter(self, X, t): pass
    
        
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

    def enter(self, X, t): pass
     
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




    
class GuidanceSearching(GuidancePurePursuit):
    def __init__(self, dm, traj, trim_args = {'h':0, 'va':12, 'gamma':0}, dt=0.01):
        GuidancePurePursuit.__init__(self, dm, p3_traj3d.CircleRefTraj(c=[0., 0., 0.], r=50.), trim_args, dt)
