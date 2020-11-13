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
        if 1:
            _zd = 1.#; zd = va*sin(gamma)
            gamma = np.arcsin(_zd/trim_args['va'])
            trim_args['gamma'] = gamma;print(gamma)
            self.Xe1, self.Ue1 = dm.trim(trim_args, report=False)
            #_thr, _ail, _ele = 0, 1, 2
            self.ff_d_thr_of_zd = -(self.Ue1[0]-self.Ue[0])
            print('ff_d_thr_of_zd {}'.format(self.ff_d_thr_of_zd))
            self.ff_d_theta_of_zd = -gamma
            print('ff_d_theta_of_zd {}'.format(self.ff_d_theta_of_zd))
        
        self.traj = traj
        self.carrot = np.zeros(3)
        self.phi_sp, self.theta_sp = self.Xe[p3_fr.SixDOFAeroEuler.sv_phi], self.Xe[p3_fr.SixDOFAeroEuler.sv_theta]
        self.att_ctl = p3_pil.AttCtl(self.Xe, self.Ue, dm, dt)
        self.T_w2b_ref = np.eye(4)
        self.sum_err_z, self.sum_err_zd, self.sum_err_v = 0, 0, 0
        self.v_mode = self.v_mode_throttle#self.v_mode_alt
        self.v_sp = trim_args['va']
        self.z_sp, self.zd_sp = 0., 0.
        self.throttle_sp = 0.

    def enter(self, X, t): pass
     
    def _compute_roll_setpoint(self, t, X, traj, max_phi=np.deg2rad(45)):
        self.b2c_ned = self.carrot-self.pos
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
        if Xee is None: print('No Xee')
        self.pos = X[p3_fr.SixDOFAeroEuler.sv_slice_pos]
        self.proj, self.carrot, self.zd_sp_traj = traj.get_proj_and_carrot(self.pos, lookahead_dist=15.)
        self._compute_roll_setpoint(t, X, traj)
        z, v = X[p3_fr.SixDOFAeroEuler.sv_z], X[p3_fr.SixDOFAeroEuler.sv_va]
        self.z_sp = self.carrot[2]
        v_sp, z_sp = self.v_sp, self.carrot[2]#Xe[p3_fr.SixDOFAeroEuler.sv_z]
        self.sum_err_v += (v-v_sp); self.sum_err_z += (z-z_sp)
        zd = Xee[p3_fr.SixDOFEuclidianEuler.sv_zd]
        self.sum_err_zd += (zd-self.zd_sp)

        # vertical control
        dthrottle, dtheta = self._v_ctl(v, z, zd, z_sp, debug=False)
        # velocity control
        self.theta_sp = self.Xe[p3_fr.SixDOFAeroEuler.sv_theta] + 2.5e-5*self.sum_err_v + dtheta
        # piloting
        U = self.att_ctl.get(t, X, self.phi_sp, self.theta_sp)
        if  self.v_mode == self.v_mode_throttle:
            U[0] = self.throttle_sp
        else:
            U[0] += dthrottle
        U[0] = np.clip(U[0], 0., 1.)             # throttle
        U[1] = np.clip(U[1], -np.pi/4, np.pi/4)  # ailerons
        U[3] = np.clip(U[3], -np.pi/4, np.pi/4)  # rudder
        return U

    def _inner_loop(self,zd, zd_sp, Kp=0.2, Ki=1e-4, Kff=0):#1e-2):
        zd_err = zd - zd_sp
        d_throttle =  Kp*zd_err + Ki*self.sum_err_zd + self.ff_d_thr_of_zd*zd_sp
        #d_throttle =  Kp*zd_err + Ki*self.sum_err_zd + Kff*zd_sp
        d_theta = self.ff_d_theta_of_zd * zd_sp
        return d_throttle, d_theta
    
    def _v_ctl(self, v, z, zd, z_sp, debug=False):
        _thr, _ail, _ele = 0, 1, 2
        # elevator compensation for banking
        #v_sp, z_sp, theta_sp = self.v_sp, self.Xe[p3_fr.SixDOFAeroEuler.sv_z], self.Xe[p3_fr.SixDOFAeroEuler.sv_theta]
        #d_ele = 0#self.sum_err_v*-0.00001 # -np.deg2rad(3.)*np.abs(self.phi_sp)/max_phi
        #U[_ele] += d_ele
        if self.v_mode == self.v_mode_throttle:
            d_throttle = self.throttle_sp
            d_theta = 0.
        elif self.v_mode == self.v_mode_vz:
            d_throttle, d_theta = self._inner_loop(zd, self.zd_sp)
        else: # v_mode_alt  # throttle integrator
            z_err = z - self.z_sp
            self.zd_sp = np.clip(-0.1*z_err, -1, 1) + np.clip(-1e-5*self.sum_err_z, -1, 1) + self.zd_sp_traj
            #d_throttle = + z_err*1e-2 # 5e-6
            #d_throttle = self.sum_err_z*0.000005
            d_throttle, d_theta = self._inner_loop(zd, self.zd_sp)
        if debug:
            fmt = 'z {:.1f} v {:.1f} m/s phi {:.1f}/{:.1f} deg theta {:.1f}/{:.1f} deg dthrottle {:.1f} % dele {:.1f} deg s_ev: {:.0f} s_ez: {:.0f}'
            print(fmt.format(z, v, np.rad2deg(phi), np.rad2deg(self.phi_sp), np.rad2deg(theta), np.rad2deg(self.theta_sp), d_throttle*100, np.rad2deg(d_ele), self.sum_err_v, self.sum_err_z))
        
        return d_throttle, d_theta

class GuidancePurePursuitLogger:
    def __init__(self):
        self.carrot = []
        self.zd_sp = []
        self.zd_sp_traj = []
        self.att_sp = []
        
    def record(self, _ctl):
        self.att_sp.append(np.array([_ctl.phi_sp, _ctl.theta_sp]))
        self.carrot.append(np.array(_ctl.carrot))
        self.zd_sp.append(_ctl.zd_sp)
        self.zd_sp_traj.append(_ctl.zd_sp_traj)

    def save(self, time, X, U, filename):
        np.savez(filename, time=time, X=X, U=U, carrot=self.carrot, zd_sp=self.zd_sp, zd_sp_traj=self.zd_sp_traj, att_sp=self.att_sp)
        print('saved {}'.format(filename))

    def load(self, filename):
        _data =  np.load(filename)
        labels = ['time', 'X', 'U', 'carrot', 'zd_sp', 'zd_sp_traj', 'att_sp']
        time, X, U, self.carrot, self.zd_sp, self.zd_sp_traj, self.att_sp = [_data[k] for k in labels]
        print('loaded {}'.format(filename))
        return time, X, U

    def plot3D(self, time, X, _ctl, atm, ref_traj=None):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        _val = np.asarray(X[:,2])#self.meas_vz)
        p3_pu.plot_3D_traj(ref_traj=ref_traj, X=X, fig=fig, ax=ax)#, val=_val)
        #p3_pu.plot_3D_wind(atm,  xspan=100, h0=0, hspan=-40, dh=-30., figure=fig, ax=ax) 
        carrot = np.array(self.carrot)
        #ax.plot(carrot[:,0], carrot[:,1], carrot[:,2], color='r', label='center')
        p3_pu.set_3D_axes_equal()


    def plot_slice_nu(self, time, X, U, _ctl, atm, ref_traj=None, n0=-40, n1=100, dn=5., h0=0., h1=60., dh=2.):
        p3_pu.plot_slice_wind_nu(atm, n0=n0, n1=n1, dn=dn, e0=0., h0=h0, h1=h1, dh=dh, zdir=-1.,
                                 show_quiver=True, show_color_bar=True, title="North Up slice",
                                 figure=None, ax=None)
        plt.plot(X[:,0],-X[:,2], color='r')
        if ref_traj is not None:
            pts = ref_traj.get_points()
            xs, ys, zs = pts[:,0], pts[:,1], pts[:,2] 
            #ax = plt.gca()
            #ax.plot(xs, ys, zs, color='g', label='ref trajectory')
            plt.plot(xs, -zs, color='g', label='reference')
        
        
        
    def plot_chronograms(self, time, X, U, _ctl, atm):
        ''' FIXME... X is Xae, is that mandatory? '''
        _s = p3_fr.SixDOFEuclidianEuler
        p1_fw_dyn.plot_trajectory_ae(time, X, U, figure=None, window_title='chronogram', legend=None, label='', filename=None, atm=atm)
        Xee = np.array([p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(_X, atm, _t) for _X, _t in zip(X, time)])

        # x,y carrot
        ax = plt.subplot(5,3,1); plt.plot(time, np.asarray(self.carrot)[:,0], label='carrot')
        ax = plt.subplot(5,3,2); plt.plot(time, np.asarray(self.carrot)[:,1], label='carrot')
        # h(up) instead of z(down)
        ax=plt.subplot(5,3,3); plt.cla()
        plt.plot(time, -Xee[:,_s.sv_z], label='plant')
        plt.plot(time, -np.asarray(self.carrot)[:,_s.sv_z], label='setpoint')
        p3_pu.decorate(ax, title="$h$", ylab="m", min_yspan=1., legend=True)
        # hdot (instead of beta)
        ax=plt.subplot(5,3,6);plt.cla()
        plt.plot(time, -Xee[:,_s.sv_zd], label="plant")
        plt.plot(time, -np.asarray(self.zd_sp), label="setpoint")
        plt.plot(time, -np.asarray(self.zd_sp_traj), label="setpoint_ff")
        p3_pu.decorate(ax, title="$\dot{h}$", ylab="m/s", min_yspan=1., legend=True)
        # attitude setpoint
        phi_sp, theta_sp = np.asarray(self.att_sp)[:,0], np.asarray(self.att_sp)[:,1] 
        ax=plt.subplot(5,3,7)
        plt.plot(time, np.rad2deg(phi_sp), label="setpoint")
        ax=plt.subplot(5,3,8)
        plt.plot(time, np.rad2deg(theta_sp), label="setpoint")

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
