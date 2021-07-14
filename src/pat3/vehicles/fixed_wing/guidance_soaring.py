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


class GuidanceAroun(p3_guid.GuidanceAuto1):
    def __init__(self, dm, traj, trim_args={'h':0, 'va':12, 'gamma':0}, dt=0.01):
        Xe, Ue = dm.trim(trim_args, report=False)
        p3_guid.GuidanceAuto1.__init__(self, Xe, Ue, dm, dt)
        self.energies, self.d_energies, self.d2_energies = [], [], []
        self.dt = dt

    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        _s = p3_fr.SixDOFAeroEuler
        alt, va = -X[_s.sv_z], X[_s.sv_va]
        self.energies.append(alt + 0.5 * va**2 / 9.81)
        d2 = 0.
        if len(self.energies) >= 2:
            self.d_energies.append((self.energies[-1]-self.energies[-2])/self.dt)
            if len(self.d_energies) >= 2:
                self.d2_energies.append((self.d_energies[-1]-self.d_energies[-2])/self.dt)
                d2 = self.d2_energies[-1]
        phi0, k = 16, 4#12
        self.phi_sp = np.deg2rad(np.clip(phi0-k*d2, -40, 40))
        self.set_setpoints(self.phi_sp, theta_sp=42.)
        return p3_guid.GuidanceAuto1.get(self, t, X, Xee, Yc, debug)

class ArounLogger():
    def __init__(self):
        self.phi_sps, self.energies, self.d_energies, self.d2_energies = [], [], [], []

    def record(self, _ctl):
        self.phi_sps.append(_ctl.phi_sp)
        self.energies.append(_ctl.energies[-1])
        try:
            self.d_energies.append(_ctl.d_energies[-1])
        except IndexError:
            self.d_energies.append(0)
        try:
            self.d2_energies.append(_ctl.d2_energies[-1])
        except IndexError:
            self.d2_energies.append(0)

    def save(self, time, X, U, filename):
        np.savez(filename, time=time, X=X, U=U, phi_sp=self.phi_sps, energy=self.energies,
                 d_energy=self.d_energies, d2_energy=self.d2_energies)
        print('saved {}'.format(filename))

    def load(self, filename):
        _data =  np.load(filename)
        labels = ['time', 'X', 'U', 'phi_sp', 'energy', 'd_energy', 'd2_energy']
        time, X, U, self.phi_sps, self.energies, self.d_energies, self.d2_energies = [_data[k] for k in labels]
        print('loaded {}'.format(filename))
        return time, X, U

    def plot_chronograms(self, time, X, _ctl, atm, title=None, window_title=None, fig=None, axes=None):
        _s = p3_fr.SixDOFEuclidianEuler
        fig = plt.figure(tight_layout=True, figsize=[6., 4.8]) if fig is None else fig
        if window_title is not None: fig.canvas.set_window_title(window_title)
        axes = fig.subplots(4, 1) if axes is None else axes
        axes[0].plot(time,  self.energies, label='energy')
        p3_pu.decorate(axes[0], title='energy', ylab='J/kg', legend=True)
        axes[1].plot(time,  self.d_energies, label='d_energy')
        p3_pu.decorate(axes[1], title='d_energy', ylab='J/kg/s', legend=True)
        axes[2].plot(time,  self.d2_energies, label='d2_energy')
        p3_pu.decorate(axes[2], title='d2_energy', ylab='J/kg/s/s', legend=True)
        axes[3].plot(time,  np.rad2deg(self.phi_sps), label='sp')
        axes[3].plot(time,  np.rad2deg(X[:,_s.sv_phi]), label='real')
        p3_pu.decorate(axes[3], title='phi', ylab='deg', legend=True)
        return fig, axes

    def plot3D(self, time, X, _ctl, atm):
        fig = plt.figure(tight_layout=True, figsize=[6., 4.8])
        ax = fig.add_subplot(111, projection='3d')
        p3_pu.plot_3D_traj(ref_traj=None, X=X, fig=fig, ax=ax)
        p3_pu.set_3D_axes_equal()
        ax.view_init(-166., -106.)
        return fig, ax

    def plot_slice_nu(self, time, X, U, _ctl, atm, n_spec, e0, h_spec,
                      title="Soaring North Up slice", window_title=None, fig=None, ax=None):
        (n0, n1, dn), (h0, h1, dh) = n_spec, h_spec
        fig, ax = p3_pu.plot_slice_wind_nu(atm, n0=n0, n1=n1, dn=dn, e0=e0, h0=h0, h1=h1, dh=dh, zdir=-1.,
                                           show_quiver=True, show_color_bar=True, title=title,
                                           figure=fig, ax=ax)
        ax.plot(X[:,0],-X[:,2])
        return fig, ax
    
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
    
    def enter(self, Xae, t): # FIXME, we need the initial circle center
        _s = p3_fr.SixDOFAeroEuler
        self.vario.reset(t=t, alt=-Xae[_s.sv_z], va=Xae[_s.sv_va], phi=Xae[_s.sv_phi])
        
    def _thermal_centering(self, t, X, Xee):
        _s = p3_fr.SixDOFAeroEuler
        X_pos = X[_s.sv_slice_pos]
        # climb rate
        self.meas_vz = Xee[p3_fr.SixDOFEuclidianEuler.sv_slice_vel][2]
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
        fig = plt.figure(tight_layout=True, figsize=[6., 4.8])
        ax = fig.add_subplot(111, projection='3d')
        _val = np.asarray(self.meas_vz)
        p3_pu.plot_3D_traj(ref_traj=None, X=X, fig=fig, ax=ax)#, val=_val)
        #p3_pu.plot_3D_wind(atm,  xspan=100, h0=0, hspan=-10, dh=-30., figure=fig, ax=ax) 
        center = np.array(self.center)
        ax.plot(center[:,0], center[:,1], center[:,2], color='r', label='center')
        p3_pu.set_3D_axes_equal()
        ax.view_init(-166., -106.)
        return fig, ax


    def plot_chronograms(self, time, X, _ctl, atm, title=None, window_title=None, fig=None, axes=None):
        _s = p3_fr.SixDOFEuclidianEuler
        fig = plt.figure(tight_layout=True, figsize=[6., 4.8]) if fig is None else fig
        if window_title is not None: fig.canvas.set_window_title(window_title)
        axes = fig.subplots(3, 1) if axes is None else axes
        center = np.array(self.center)

        axes[0].plot(time,  center[:,0], label='est')
        axes[1].plot(time,  center[:,1], label='est')

        try:   # plot center of atmosphere for analytic model that have them
            foo = atm.radius
            axes[0].plot(time,  np.ones(len(time))*atm.center[0], label='truth')
            axes[1].plot(time,  np.ones(len(time))*atm.center[1], label='truth')
        except AttributeError: pass
        p3_pu.decorate(axes[1], title='cy', ylab='m', legend=True)
        p3_pu.decorate(axes[0], title='cx', ylab='m', legend=True)

        axes[2].plot(time, self.meas_netto, label='netto')
        axes[2].plot(time, self.meas_vz, label='vz')
        p3_pu.decorate(axes[2], title='measurement (vz)', ylab='m/s', legend=True)
        return fig, axes

    def plot_slice_nu(self, time, X, U, _ctl, atm, n_spec, e0, h_spec,
                      title="Soaring North Up slice", window_title=None, fig=None, ax=None):
        (n0, n1, dn), (h0, h1, dh) = n_spec, h_spec
        fig, ax = p3_pu.plot_slice_wind_nu(atm, n0=n0, n1=n1, dn=dn, e0=e0, h0=h0, h1=h1, dh=dh, zdir=-1.,
                                           show_quiver=True, show_color_bar=True, title=title,
                                           figure=fig, ax=ax)
        ax.plot(X[:,0],-X[:,2])
        center_ned = np.array(self.center)
        ax.plot(center_ned[:,0],-center_ned[:,2], 'r')
        return fig, ax
