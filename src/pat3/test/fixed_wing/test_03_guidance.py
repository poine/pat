#! /usr/bin/env python
import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.atmosphere as p3_atm
import pat3.vehicles.rotorcraft.multirotor_trajectory as p3_mt
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.algebra as pal

import control.matlab

import test_02_att_ctl


class RefTraj:
    def __init__(self):
        pass

    def points(self):
        return None
    

class CircleRefTraj:
    def __init__(self, c=(0, 0, 10), r=40.):
        self.c, self.r = np.asarray(c), r

        alphas = np.linspace(0, 2*np.pi, 360)
        xs = c[0] + self.r*np.cos(alphas)
        ys = c[1] + self.r*np.sin(alphas)
        zs = c[2] + 5*np.sin(alphas)#np.zeros(len(alphas))
        self.points = np.vstack((xs, ys, zs)).T
        self.dists = alphas * r

    def get_alpha(self, p):
        cp = p-self.c
        return np.arctan2(cp[1], cp[0])
        
        
    def get_points(self):
        return self.points

    def get_closest_point(self, p):
        cp = p-self.c
        return self.c + cp/np.linalg.norm(cp)

    def get_point_ahead(self, p, lookahead_dist):
        alpha0 = self.get_alpha(p)
        alpha1 = alpha0 + lookahead_dist/self.r
        p1 = self.c+[self.r*np.cos(alpha1), self.r*np.sin(alpha1), 5*np.sin(alpha1)]
        return p1


class CircleRefTraj2(p3_mt.Circle):
    def __init__(self, c=(0, 0, 10), r=40.):
        p3_mt.Circle.__init__(self, c=[0, 0, 0], r=r, v=15., alpha0=0, dalpha=2*np.pi)

    #def _zfun(self, alpha): return  5*np.sin(alpha) 
    def _zfun(self, alpha): return  np.zeros_like(alpha)
        
    def get_points(self):
        alphas = np.linspace(self.alpha0, self.alpha0+self.dalpha, 360*self.dalpha/np.pi/2)
        xs = self.c[0] + self.r*np.cos(alphas)
        ys = self.c[1] + self.r*np.sin(alphas)
        zs = self.c[2] + self._zfun(alphas)#np.zeros(len(alphas))
        return np.vstack((xs, ys, zs)).T

        
    def get_alpha(self, p): # FIXME alpha0
        cp = p-self.c
        return np.arctan2(cp[1], cp[0])
    
    def get_point_ahead(self, p0, lookahead_dist):
        alpha0 = self.get_alpha(p0)
        alpha1 = alpha0 + lookahead_dist/self.r
        p1 = self.c+[self.r*np.cos(alpha1), self.r*np.sin(alpha1), self._zfun(alpha1)]
        #p0p1 = p1-p0
        #print('alpha0 {:.3f} alpha1 {:.3f} p0 {} p1 {} dist {:.2f}'.format(alpha0, alpha1, p0, p1, np.linalg.norm(p0p1)))
        return p1

    
    
class SpiralRefTraj:
    def __init__(self, c=(0, 0, 10), r=30.):
        pass

def norm_mpi_pi(v): return ( v + np.pi) % (2 * np.pi ) - np.pi

class Guidance:
    def __init__(self, dm, traj, trim_args = {'h':0, 'va':12, 'gamma':0}, dt=0.01):
        self.Xe, self.Ue = dm.trim(trim_args, debug=True)
        self.traj = traj
        self.att_ctl = test_02_att_ctl.AttCtl(self.Xe, self.Ue, dm, dt)
        self.spv_size = 2
        self.T_w2b_ref = np.eye(4)
        self.sum_err_z, self.sum_err_v = 0, 0

    def get(self, t, X, Yc=None, debug=False):
        my_pos = X[p1_fw_dyn.sv_slice_pos]
        self.carrot = self.traj.get_point_ahead(my_pos, 15.)
        self.b2c_ned = self.carrot-my_pos
        R_w2b = pal.rmat_of_euler(X[p1_fw_dyn.sv_slice_eul])
        self.b2c_b = np.dot(R_w2b, self.b2c_ned)
        z, v = X[p1_fw_dyn.sv_z], X[p1_fw_dyn.sv_v]
        phi, theta = X[p1_fw_dyn.sv_phi], X[p1_fw_dyn.sv_theta]
        if 0:
            self.R = (np.linalg.norm(self.b2c_b)**2)/(2*self.b2c_b[1])
            # R.g.tan(phi) = v^2
            phi_sp = np.arctan(v**2/self.R/9.81)
        else:
            self.R = 0
            err_psi = X[p1_fw_dyn.sv_psi] - np.arctan2(self.b2c_ned[1], self.b2c_ned[0])
            err_psi = norm_mpi_pi(err_psi)
            phi_sp = -0.75*err_psi

        max_phi = np.deg2rad(45)
        self.phi_sp = np.clip(phi_sp, -max_phi, max_phi)
        self.theta_sp = self.Xe[p1_fw_dyn.sv_theta] + 0.000025*self.sum_err_v
        U = self.att_ctl.get(t, X, self.phi_sp, self.theta_sp)
        _thr, _ail, _ele = 0, 1, 2
        # elevator compensation for banking
        v_sp, z_sp, theta_sp = self.Xe[p1_fw_dyn.sv_v], self.Xe[p1_fw_dyn.sv_z], self.Xe[p1_fw_dyn.sv_theta]
        self.sum_err_v += (v-v_sp)
        self.sum_err_z += (z-z_sp)
        d_ele = 0#self.sum_err_v*-0.00001 # -np.deg2rad(3.)*np.abs(phi_sp)/max_phi
        U[_ele] += d_ele
        # throttle integrator
        d_throttle = self.sum_err_z*0.000005
        U[_thr] += d_throttle
        if debug:
            fmt = 'z {:.1f} v {:.1f} m/s phi {:.1f}/{:.1f} deg theta {:.1f}/{:.1f} deg dthrottle {:.1f} % dele {:.1f} deg s_ev: {:.0f} s_ez: {:.0f}'
            print(fmt.format(z, v, np.rad2deg(phi), np.rad2deg(self.phi_sp), np.rad2deg(theta), np.rad2deg(self.theta_sp), d_throttle*100, np.rad2deg(d_ele), self.sum_err_v, self.sum_err_z))
        
        return U#self.Ue

def run_simulation(dm, ref_traj, tf=60.5, dt=0.01, trim_args = {'h':0, 'va':12, 'gamma':0}, plot=False):
    time = np.arange(0, tf, dt)
    X = np.zeros((len(time), dm.sv_size))
    U = np.zeros((len(time),  dm.input_nb()))
    carrots, att_sp = np.zeros((len(time),  3)), np.zeros((len(time), 2))
    ctl = Guidance(dm, ref_traj, trim_args, dt)
    X[0] = dm.reset(ctl.Xe)
    for i in range(1, len(time)):
        U[i-1] = ctl.get(time[i-1], X[i-1])
        carrots[i-1] = ctl.carrot; att_sp[i-1] = ctl.phi_sp, ctl.theta_sp 
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1])
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

    
def plot_3d_traj(ref_traj, X=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    pts = ref_traj.get_points()
    xs, ys, zs = pts[:,0], pts[:,1], pts[:,2] 
    ax.plot(xs, ys, zs, color='g', label='ref trajectory')
    ax.plot(xs, ys, np.zeros(len(xs)), linestyle='--', color='k', linewidth=0.5, label='ground track')
    # https://stackoverflow.com/questions/36013063/what-is-the-purpose-of-meshgrid-in-python-numpy
    if 1:
        verts = []
        for i in range(1, len(pts)):
            verts.append([(xs[i-1], ys[i-1], 0), (xs[i-1], ys[i-1], zs[i-1]),
                          (xs[i], ys[i], zs[i]), (xs[i], ys[i], 0)])
        mesh = Poly3DCollection(verts, alpha=0.2, linewidth=0, facecolor=[0.5, 0.5, 1])
        ax.add_collection3d(mesh)

    if X is not None:
        ax.plot(X[:,0], X[:,1], X[:,2], color='b', label='aircraft trajectory')
        
    p3_pu.set_3D_axes_equal()
    plt.legend(loc='best')

def plot_3d_wind(atm):
    ax = plt.gca()
    # Make the grid
    x, y, z = np.meshgrid(np.linspace(-30., 30., 20),
                          np.linspace(-30., 30., 20),
                          np.linspace(0., 10, 10))

    # Make the direction data for the arrows
    u = np.sin(np.pi * x) * np.cos(np.pi * y) * np.cos(np.pi * z)
    v = -np.cos(np.pi * x) * np.sin(np.pi * y) * np.cos(np.pi * z)
    w = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * x) * np.cos(np.pi * y) *
         np.sin(np.pi * z))
    
    ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True)

    
    
def main(param_filename, trim_args = {'h':0, 'va':11, 'gamma':0}):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    ref_traj = CircleRefTraj2(c=[10, 0, 10])
    #time, X, U = test_02_pitch_ctl.run_simulation(dm, plot=False)
    time, X, U = run_simulation(dm, ref_traj, plot=True)
    
    #atm =  p3_atm.Atmosphere()
    plot_3d_traj(ref_traj, X)
    #plot_3d_wind(atm)
    plt.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    dirname, filename = os.path.split(os.path.abspath(__file__))
    param_filename = os.path.abspath(os.path.join(dirname, '../../../../data/vehicles/cularis.xml'))
    main(param_filename)
