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

import pat3.vehicles.fixed_wing.guidance as p3_guid
import pat3.vehicles.fixed_wing.guidance_ardusoaring as p3_guidardu
    
class GuidanceSoaring(p3_guid.GuidancePurePursuit):

    def __init__(self, dm, traj, trim_args={'h':0, 'va':12, 'gamma':0}, dt=0.01):
        p3_guid.GuidancePurePursuit.__init__(self, dm, p3_traj3d.CircleRefTraj(c=[0., 0., 0.], r=15.), trim_args, dt)
        self.v_mode_throttle = p3_guid.GuidancePurePursuit.v_mode_throttle
        self.throttle_sp = 0.
        self.k_centering = 0.002
        # circle discretization
        self.nb_angles = 360
        self.avg_climb = np.zeros(self.nb_angles)
        self.avg_climb[:] = np.float('nan')
        self.meas_vz = 0.
        
        self.meas_netto = 0.
        self.vario = p3_guidardu.NettoVario()

    def _idx_of_angle(self, angle): return int(p3_u.norm_0_2pi(angle)/2/np.pi*float(self.nb_angles))
    def _angle_of_idx(self, idx): return p3_u.norm_mpi_pi(idx/float(self.nb_angles)*2*np.pi)    

    def set_param(self, _v): self.k_centering = _v
    def set_radius(self, _v): self.traj.r = _v
    def set_circle_center(self, _c): self.traj.c = np.asarray(_c)
    
    def enter(self, Xae, t):
        _s = p3_fr.SixDOFAeroEuler
        alt, va, phi, psi = -Xae[_s.sv_z], Xae[_s.sv_va], Xae[_s.sv_phi], Xae[_s.sv_psi]
        self.vario.reset(t=t, alt=alt, va=va, phi=phi)
        
    def _thermal_centering(self, t, X, Xee):
        _s = p3_fr.SixDOFAeroEuler
        X_pos = X[_s.sv_slice_pos]
        # climb rate
        self.meas_vz = Xee[p3_fr.SixDOFEuclidianEuler.sv_slice_vel]
        # specific energy rate
        alt, va, phi = -X[_s.sv_z], X[_s.sv_va], X[_s.sv_phi] 
        self.vario.update(t, alt, va, phi)
        self.meas_netto = -self.vario.reading
        alpha = self.traj.get_alpha(X_pos)
        idx_alpha = self._idx_of_angle(alpha)
        lp = 0.25
        self.avg_climb[idx_alpha] = self.meas_netto#self.meas_vz#lp*self.avg_climb[idx_alpha] + (1-lp)*ivel_ned[2]
        if np.isnan(self.avg_climb).any():
            self.delta_center = np.zeros(3)
        else:
            amin, amax = np.argmin(self.avg_climb), np.argmax(self.avg_climb)
            alpha_min, alpha_max = self._angle_of_idx(amin), self._angle_of_idx(amax)
            d_climb = self.avg_climb[amax] - self.avg_climb[amin]
            alpha_max_grad = alpha_max#(alpha_max+alpha_min + np.pi)/2
            self.delta_center = np.array([np.cos(alpha_max_grad), np.sin(alpha_max_grad), 0.])*d_climb*self.k_centering
            
    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        # offset circle toward vz gradient
        self._thermal_centering(t, X, Xee)
        self.traj.offset(self.delta_center)
        
        self.traj.c[2] = Xee[p3_fr.SixDOFEuclidianEuler.sv_z]
        return p3_guid.GuidancePurePursuit.get(self, t, X, Xee, Yc, debug)




class Logger:
    def __init__(self):
        self.center = []
        self.meas_vz = []
        self.meas_netto = []
        
    def record(self, _ctl):
        self.center.append(np.array(_ctl.traj.c))
        self.meas_vz.append(_ctl.meas_vz)
        self.meas_netto.append(_ctl.meas_netto)
    
    def save(self, time, X, U, filename):
        np.savez(filename, time=time, X=X, U=U, center=self.center, meas_vz=self.meas_vz, meas_netto=self.meas_netto)
        print('saved {}'.format(filename))

    def load(self, filename):
        _data =  np.load(filename)
        labels = ['time', 'X', 'U', 'center', 'meas_vz', 'meas_netto']
        time, X, U, self.center, self.meas_vz, self.meas_netto = [_data[k] for k in labels]
        print('loaded {}'.format(filename))
        return time, X, U

    
    def plot3D(self, time, X, _ctl, atm):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        _val = np.asarray(meas_vz)
        p3_pu.plot_3D_traj(ref_traj=None, X=X, fig=fig, ax=ax)#, val=_val)
        p3_pu.plot_3D_wind(atm,  xspan=100, h0=0, hspan=-10, dh=-30., figure=fig, ax=ax) 
        center = np.array(self.center)
        ax.plot(center[:,0], center[:,1], center[:,2], color='r', label='center')
        p3_pu.set_3D_axes_equal()


    def plot_chronograms(self, time, X, _ctl, atm):
        _s = p3_fr.SixDOFEuclidianEuler
        plt.figure()
        center = np.array(self.center)
        meas_vz, meas_netto = np.array(self.meas_vz), np.array(self.meas_netto)
        plots = [("cx", "m", center[:,0]),
                 ("cy", "m", center[:,1])]
        for i, (title, ylab, data) in enumerate(plots):
            ax = plt.subplot(3,1,i+1)
            plt.plot(time, data)
            p3_pu.decorate(ax, title=title, ylab=ylab)#, legend=True)
        ax = plt.subplot(3,1,1); plt.plot(time,  np.ones(len(time))*atm.center[0], label='truth')
        ax = plt.subplot(3,1,2); plt.plot(time,  np.ones(len(time))*atm.center[1], label='truth')
        ax = plt.subplot(3,1,3)
        plt.plot(time,  meas_netto, label='netto')
        plt.plot(time, meas_vz, label='vz')
        p3_pu.decorate(ax, title='measurement', ylab='m/s', legend=True)
