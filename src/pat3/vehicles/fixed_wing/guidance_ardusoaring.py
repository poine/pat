import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt
import matplotlib
import mpl_toolkits

'''
This a test re implementation of ardusoaring
'''


import pat3.vehicles.fixed_wing.guidance as p3_guid
import pat3.trajectory_3D as p3_traj3d
import pat3.frames as p3_fr
import pat3.atmosphere as p3_atm
import pat3.plot_utils as p3_pu

class NettoVario:
    VA_FLT_K = 0.05
    TE_FLT_K = 0.03
    DISP_TE_FLT_K = 0.15
    GRAVITY_MSS = 9.80665
    
    def __init__(self):
        # values obtained from fit_netto_vario.py 
        self.polar_K   =  49.5050092764    # 25.6   # Cl factor 2*m*g/(rho*S) @Units: m.m/s/s
        self.polar_CD0 =   0.0122440667444 #  0.027 # Zero lift drag coef
        self.polar_B   =   0.0192172535765 #  0.031 # Induced drag coeffient
        self.reset()

    def reset(self, t=0, alt=0, va=10, phi=0):
        self.va_flt = va
        self.reading, self.reading_flt, self.reading_disp = 0., 0., 0.
        self._store_state(t, alt, va, phi, alt + 0.5 * va**2 / self.GRAVITY_MSS)
        self.reading = self.reading_flt = self.reading_disp = self.correct_netto_rate(0.0, phi, va)
        
        
    def _store_state(self, t, alt, va, phi, total_E): 
        self.last_t = t
        self.last_alt = alt
        self.last_va = va
        self.last_phi = phi
        self.last_total_E = total_E
        
    def update(self, t, alt, va, phi):
        #if np.abs(alt-self.last_alt):
        self.va_flt = self.VA_FLT_K * va + (1.-self.VA_FLT_K) * self.va_flt    # low pass airspeed
        total_E = alt + 0.5 * self.va_flt**2 / self.GRAVITY_MSS                # compute total energy
        sinkrate = self.correct_netto_rate(0.0, (phi + self.last_phi) / 2., self.va_flt)
        d_E, d_t = total_E - self.last_total_E, t-self.last_t
        if d_t>1e-12:
            self.reading = d_E/d_t + sinkrate                                      # unfiltered netto rate
            self.reading_flt = self.TE_FLT_K * self.reading + (1.-self.TE_FLT_K) * self.reading_flt
            self.reading_disp = self.DISP_TE_FLT_K * self.reading + (1.-self.DISP_TE_FLT_K) * self.reading_disp
            self._store_state(t, alt, va, phi, total_E)
            
    def correct_netto_rate(self, climb_rate, phi, va):
        CL0 =  self.polar_K / va**2
        C1 = self.polar_CD0 / CL0  # constant describing expected angle to overcome zero-lift drag
        C2 = self.polar_B * CL0    # constant describing expected angle to overcome lift induced drag at zero bank               
        #cosphi = (1 - phi * phi / 2)                                  # first two terms of mclaurin series for cos(phi)
        netto_rate = climb_rate + va * (C1 + C2 / np.cos(phi)**2)#(cosphi * cosphi))  # effect of aircraft drag removed
        # Remove acceleration effect - needs to be tested.
        # float temp_netto = netto_rate;
        # float dVdt = SpdHgt_Controller->get_VXdot();
        # netto_rate = netto_rate + aspd*dVdt/GRAVITY_MSS;
        # gcs().send_text(MAV_SEVERITY_INFO, "%f %f %f %f",temp_netto,dVdt,netto_rate,barometer.get_altitude());
        return netto_rate

class ThermalFilter:
    sv_strength = 0
    sv_radius   = 1
    sv_x        = 2
    sv_y        = 3
    sv_size     = 4

    _q_strength = 0.001
    _q_pos      = 0.03

    _init_thermal_dist = 30.#5.
    INITIAL_THERMAL_STRENGTH = 2.0
    INITIAL_THERMAL_RADIUS = 30.0
    INITIAL_STRENGTH_COVARIANCE = 0.0049
    INITIAL_RADIUS_COVARIANCE = 50**2     # 2500
    INITIAL_POSITION_COVARIANCE = 18**2   #  324
    
    def __init__(self):
        self.q1, self.q2 =  0.001**2, 0.03**2
        self.Q = np.diag((self.q1, self.q2, self.q2, self.q2))
        self.R = 0.45**2
        self.reset(yaw=0.)
        
    def reset(self, yaw, xc0=None, yc0=None, s0=INITIAL_THERMAL_STRENGTH, r0=INITIAL_THERMAL_RADIUS):
        xc0 = xc0 if xc0 is not None else self._init_thermal_dist * np.cos(yaw)
        yc0 = yc0 if yc0 is not None else self._init_thermal_dist * np.sin(yaw)
        self.X = np.array([s0, r0, xc0, yc0])
        self.P = np.diag([self.INITIAL_STRENGTH_COVARIANCE, self.INITIAL_RADIUS_COVARIANCE,
                          self.INITIAL_POSITION_COVARIANCE, self.INITIAL_POSITION_COVARIANCE])
        self.inov = 0
        print('ekf reset {}'.format(self.X))


    def predict_measurement(self):
        _s, _r, _x, _y = self.X
        expon = np.exp(-(_x**2 + _y**2) / _r**2)
        w = _s*expon
        H = np.zeros((1,4))
        H[0,0] = expon
        H[0,1] =  2.*_s*((_x**2+_y**2)/_r**3)*expon
        H[0,2] = -2.*(_s*_x/_r**2)*expon
        #H[0,3] = H[0,2]*_y/_x # dangerous if x is zero
        H[0,3] = -2.*(_s*_y/_r**2)*expon
        return w, H
        
    def update(self, t, z, vx, vy, disable_correction=False):
        if 1: # Prediction
            self.X[self.sv_x] -= vx
            self.X[self.sv_y] -= vy
            # predict covariance
            # P = A*ekf.P*A'+ekf.Q;
            self.P += self.Q
        # Correction
        w, H = self.predict_measurement()
        self.inov = z-w
        if not disable_correction: 
            PHt = np.dot(self.P, H.T)
            S = np.dot(H, PHt) + self.R
            K = PHt/S
            self.X += (K*self.inov).squeeze()
            #pdb.set_trace()
            #self.P = (1-np.dot(K, H))*self.P
            self.P = (np.eye(4)-np.dot(K, H))*self.P
            # P.force_symmetry()
            #print('{} ekf:  meas({}/{})'.format(0, z, w))
        else:
            print('{} ekf: correction disabled meas({}/{})'.format(0, z, w))
            #print('ekf update {}'.format(self.X))

    def center_ned(self, Xac, Xflt=None):
        Xflt = Xflt if Xflt is not None else self.X
        return Xac[p3_fr.SixDOFEuclidianEuler.sv_slice_pos] + [Xflt[self.sv_x], Xflt[self.sv_y], 0.]
        
        
class GuidanceArduSoaring(p3_guid.GuidancePurePursuit):

    def __init__(self, dm, traj, trim_args = {'h':0, 'va':12, 'gamma':0}, dt=0.01):
        self.r = 20
        p3_guid.GuidancePurePursuit.__init__(self, dm, p3_traj3d.CircleRefTraj(c=[0., 15., 0.], r=self.r), trim_args, dt)
        self.vario = NettoVario()
        self.ekf = ThermalFilter()
        
    def _store_state(self, X, t):
        self.prev_state = np.array(X)
        self.prev_t = t
        
    def enter(self, X, t):
        self._store_state(X, t)
        _s = p3_fr.SixDOFAeroEuler
        alt, va, phi, psi = -X[_s.sv_z], X[_s.sv_va], X[_s.sv_phi], X[_s.sv_psi]
        self.vario.reset(t=t, alt=alt, va=va, phi=phi)
        x0, y0 = None, None#-20., 5.
        s0, r0 = -1.5, 40
        self.ekf.reset(psi+np.pi, x0, y0, s0, r0)
        print('ekf reset: center_ned {}'.format(self.ekf.center_ned(X)))

    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        _s = p3_fr.SixDOFAeroEuler
        alt, va, phi = -X[_s.sv_z], X[_s.sv_va], X[_s.sv_phi] 
        self.vario.update(t, alt, va, phi)
        est_wind_vne = np.array([0, 0])
        # get wind corrected drift
        _s = p3_fr.SixDOFAeroEuler
        d_ne = X[_s.sv_x:_s.sv_y+1] - self.prev_state[_s.sv_x:_s.sv_y+1]
        d_t = t - self.prev_t
        d_wind_ne =  est_wind_vne * d_t
        d_ne_corrected = d_ne - d_wind_ne
        self._store_state(X, t)
        # FIXME sign baro...
        if t>3.:
            self.ekf.update(t, -self.vario.reading, d_ne_corrected[0], d_ne_corrected[1], disable_correction=False)#t<3.) # netto baro is bad at first
        #print('arduguidance update: c_ned {}'.format(self.ekf.center_ned(X)))

        if 1: # close the loop
            #t0 = 60
            #if t > t0: self.traj.r = self.r*(1+0.5*np.sin(0.1*(t-t0)))
            #dx = 10*np.sin(0.05*(t-t0)) if t > t0 else 0
            #dx = 0
            #self.traj.set_center([self.ekf.X[ThermalFilter.sv_x]+dx, self.ekf.X[ThermalFilter.sv_y], -alt])
            self.traj.set_center(self.ekf.center_ned(Xee))
        return p3_guid.GuidancePurePursuit.get(self, t, X, Xee, Yc, debug)




class Logger:
    def __init__(self):
        self.Xs, self.Ps, self.inovs = [], [], []
        self.vs_reading, self.vs_flt, self.vs_disp = [], [], []

    def record(self, _ctl):
        self.Xs.append(np.array(_ctl.ekf.X))
        self.Ps.append(_ctl.ekf.P)
        self.inovs.append(_ctl.ekf.inov)
        self.vs_reading.append(_ctl.vario.reading)
        self.vs_flt.append(_ctl.vario.reading_flt)
        self.vs_disp.append(_ctl.vario.reading_disp)

    def save(self, time, X, U, filename):
        np.savez(filename, time=time, X=X, U=U, ekfX=self.Xs, ekfP=self.Ps, inovs=self.inovs, vs_reading=self.vs_reading, vs_flt=self.vs_flt, vs_disp=self.vs_disp)
        print('saved {}'.format(filename))

    def load(self, filename):
        _data =  np.load(filename)
        labels = ['time', 'X', 'U', 'ekfX', 'ekfP', 'inovs', 'vs_reading', 'vs_flt', 'vs_disp']
        time, X, U, self.Xs, self.Ps, self.inovs, self.vs_reading, self.vs_flt, self.vs_disp = [_data[k] for k in labels]
        print('loaded {}'.format(filename))
        return time, X, U
    
        
    def plot3D(self, time, X, _ctl, atm):

        # 3D
        #p3_pu.plot_3D_traj(ref_traj=None, X=X)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('North');ax.set_ylabel('East');ax.set_zlabel('Down')        
        p3_pu.plot_3D_traj(ref_traj=None, X=X, fig=fig, ax=ax)#, val=_val)
        
        #ax.plot(center[:,0], center[:,1], center[:,2], color='r', label='center')
        # ax, fig = plt.gca(), plt.gcf()
        # ac_ned = X[:,p3_fr.SixDOFEuclidianEuler.sv_slice_pos]
        # ac_enu = p3_fr.ned_to_enu(ac_ned).reshape(-1, 1, 3)
        # _ac_ned = ac_ned.reshape(-1, 1, 3)
        # segments_ned = np.concatenate([_ac_ned[:-1], _ac_ned[1:]], axis=1)
        
        # _val = np.asarray(self.vs_flt[:-1])
        # #norm = plt.Normalize(_val.min(), _val.max())
        # norm = plt.Normalize(0., 1.5)
        
        # lc = mpl_toolkits.mplot3d.art3d.Line3DCollection(segments_ned, cmap=plt.get_cmap('viridis'), norm=norm)
        # lc.set_array(_val) 
        # lc.set_linewidth(2)
        # #pdb.set_trace()
        # ax.add_collection3d(lc)#, zs=z, zdir='z')
        # p3_pu.set_3D_axes_equal()
        
        # if 1: # center
        c_ned = np.array([_ctl.ekf.center_ned(Xac, Xflt) for Xac, Xflt in zip(X, self.Xs)])
        #     #c_ned = np.array([_ctl.ekf.center_ned(_X) for _X in X])
        ax.plot(c_ned[:,0], c_ned[:,1], c_ned[:,2], color='r', label='center')
        # if 1: #
        #     pass
        # plt.legend()

        #atm = p3_atm.AtmosphereWharington()
        #atm = p3_atm.AtmosphereThermal()
        #atm = p3_atm.AtmosphereThermal1()
        #p3_pu.plot_3D_wind(atm,  xspan=100, h0=0, hspan=10, dh=10., figure=fig, ax=ax)
        ax.view_init(-166., -106.)

    def plot2D(self, time, X, _ctl, atm, n_spec=(-100, 100, 10), e0=0, h_spec=(180., 320., 10.)):
        # 2D
        (n0, n1, dn), (h0, h1, dh) = n_spec, h_spec
        fig, ax = p3_pu.plot_slice_wind_nu(atm, n0, n1, dn, e0, h0, h1, dh)
        #, n0=-80, n1=60, dn=10., h0=180., h1=320., dh=10.
        ac_ned = X[:,p3_fr.SixDOFEuclidianEuler.sv_slice_pos]
        ac_enu = p3_fr.ned_to_enu(ac_ned).reshape(-1, 1, 3)
        #plt.plot(ac_ned[:,0], -ac_ned[:,2])
        segments_enu = np.concatenate([ac_enu[:-1], ac_enu[1:]], axis=1)

        _val = np.asarray(self.vs_flt[:-1])
        #norm = plt.Normalize(_val.min(), _val.max())
        norm = plt.Normalize(0., 1.6)

        lc = matplotlib.collections.LineCollection(segments_enu[:,:,1:], cmap='viridis', norm=norm) # north and up
        lc.set_array(_val)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        cbar = fig.colorbar(line, ax=ax)
        cbar.ax.set_ylabel('m/s (up)', rotation=270); cbar.ax.set_xlabel('netto vario')

        c_ned = np.array([_ctl.ekf.center_ned(Xac, Xflt) for Xac, Xflt in zip(X, self.Xs)])
        plt.plot(c_ned[:,0], -c_ned[:,2], 'r') # north up
        
        #pdb.set_trace()
        
    def plot_chronograms(self, time, X, _ctl, atm):
        _s = p3_fr.SixDOFEuclidianEuler
        # compute truth
        try:
            X_flt_truth = np.zeros((len(time), ThermalFilter.sv_size))
            X_flt_truth[:, ThermalFilter.sv_strength] = atm.strength
            X_flt_truth[:, ThermalFilter.sv_radius] = atm.radius
            for i in range(len(time)):
                X_flt_truth[i, ThermalFilter.sv_x:ThermalFilter.sv_y+1] = (atm.center - X[i,_s.sv_slice_pos])[:2] # FIXME... time?
        except AttributeError:
            pass
        plt.figure()
        ax = plt.subplot(2,1,1)
        plt.plot(time, self.vs_reading, label='vario')
        plt.plot(time, self.vs_flt,     label='vario_flt')
        plt.plot(time, self.vs_disp,    label='vario_disp')
        p3_pu.decorate(ax, title='netto vario', ylab='m/s', legend=True)
        ax = plt.subplot(2,1,2)
        plt.plot(time, self.inovs, label='innovation')
        p3_pu.decorate(ax, title='innovation', ylab='m/s', legend=True)
        plt.figure()
        Xs, Ps = np.array(self.Xs), np.array(self.Ps)
        plots = [("strength", "m/s", Xs[:,ThermalFilter.sv_strength]),
                 ("P_s", "m/s**2",   Ps[:,ThermalFilter.sv_strength,ThermalFilter.sv_strength]),
                 ("radius",   "m",   Xs[:,ThermalFilter.sv_radius]),
                 ("P_r", "m**2",     Ps[:,ThermalFilter.sv_radius,ThermalFilter.sv_radius]),
                 ("xc",   "m",       Xs[:,ThermalFilter.sv_x]),
                 ("P_x", "m**2",     Ps[:,ThermalFilter.sv_x,ThermalFilter.sv_x]),
                 ("yc",   "m",       Xs[:,ThermalFilter.sv_y]),
                 ("P_y", "m**2",     Ps[:,ThermalFilter.sv_y,ThermalFilter.sv_y])]
        for i, (title, ylab, data) in enumerate(plots):
            ax = plt.subplot(4,2,i+1)
            plt.plot(time, data)
            p3_pu.decorate(ax, title=title, ylab=ylab)#, legend=True)
        ax = plt.subplot(4,2,1); plt.plot(time,  X_flt_truth[:,ThermalFilter.sv_strength])
        ax = plt.subplot(4,2,3); plt.plot(time,  X_flt_truth[:,ThermalFilter.sv_radius])
        ax = plt.subplot(4,2,5); plt.plot(time,  X_flt_truth[:,ThermalFilter.sv_x])
        ax = plt.subplot(4,2,7); plt.plot(time,  X_flt_truth[:,ThermalFilter.sv_y])
        

    
