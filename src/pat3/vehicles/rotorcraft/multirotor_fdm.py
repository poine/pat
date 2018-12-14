#! /usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

import pat3.algebra as pal
import pat3.plot_utils as ppu
import pdb

'''
State vector components
'''
sv_x    = 0  # position in euclidian world
sv_y    = 1  #
sv_z    = 2  # velocity in euclidian world
sv_xd   = 3  #
sv_yd   = 4  #
sv_zd   = 5  #
sv_qi   = 6  # world to body rotation quaternion
sv_qx   = 7  #
sv_qy   = 8  #
sv_qz   = 9  #
sv_p   = 10  # rotational velocities in body frame
sv_q   = 11  # ( front, right, down )
sv_r   = 12  #
sv_size = 13

sv_slice_pos  = slice(sv_x, sv_z+1)
sv_slice_vel  = slice(sv_xd, sv_zd+1)
sv_slice_quat = slice(sv_qi, sv_qz+1)
sv_slice_rvel = slice(sv_p,sv_r+1)



'''
input vector components
'''
iv_fr   = 0
iv_br   = 1
iv_bl   = 2
iv_fl   = 3
iv_size = 4


class FDM:
    ''' An object encapsulating the dynamic model of the multirotor '''
    def __init__(self, dt=0.005):
        self.dt = dt
        self.P = Param()
        self.t, self.X = 0., np.zeros((sv_size))
        # byproducts
        self.T_w2b = np.eye(4)  # world (ned) to body (frd) homogeneous transform
        
    def trim(self):
        return trim(self.P)
    
    def reset(self, X0, t0):
        self.X, self.t = X0, t0
        self.update_byproducts()
        return self.X
    
    def run(self, tf, U):
        remaining_to_tf = tf - self.t
        while remaining_to_tf > 0:
            dt = min(self.dt, remaining_to_tf)
            #print(' integrate fdm from {} to {}'.format(self.t, self.t+dt))
            self.X = disc_dyn(self.X, self.t, U, dt, self.P)
            self.t += dt
            remaining_to_tf = tf - self.t
        self.update_byproducts()
        return self.X
        
    def update_byproducts(self):
        self.T_w2b[:3,:3] = pal.rmat_of_quat(self.X[sv_slice_quat]).T # that is freaking weird....
        self.T_w2b[:3,3] = self.X[sv_slice_pos]
        
    def plot(self, time, X, U=None, figure=None, window_title="Trajectory"):
        return plot(time, X, U, figure, window_title)


class Param():
    ''' Represents the configuration of the multirotor (mass, inertia, rotors location etc...) '''
    def __init__(self):
        self.m = 0.5
        self.J = np.diag([0.01, 0.01, 0.05])
        self.g = 9.81
        # position of rotors in body frame
        self.rotor_pos = ([0.1, 0.1, 0], [-0.1, 0.1, 0], [-0.1, -0.1, 0], [0.1, -0.1, 0])
        self.l = 0.14
        # direction of rotors ( 1 cw, -1 ccw )
        self.rotor_dir = ( -1., 1., -1., 1. )
        # torque over thrust coefficient
        self.k = 0.1
        # drag - needs love
        self.Cd = 0.2
        # precompute
        self.invJ = np.linalg.inv(self.J)
        




def trim(P):
    Xe = np.zeros(sv_size); Xe[sv_qi] = 1.
    Ue = np.ones(iv_size)*P.m*P.g / iv_size
    return Xe, Ue



def solid_cont_dyn(X, F_b, M_b, P):
    Xd = np.zeros(sv_size)
    p_w, v_w, q_w2b, om_b = X[sv_slice_pos], X[sv_slice_vel], X[sv_slice_quat], X[sv_slice_rvel]
    # Translational kinematics
    Xd[sv_slice_pos] = v_w
    # Newton for forces
    R_w2b =  pal.rmat_of_quat(q_w2b)
    Xd[sv_slice_vel] = 1./P.m*(np.dot(R_w2b.T, F_b) + [0, 0, P.m*P.g])
    # Rotational kinematics
    Xd[sv_slice_quat] = pal.quat_derivative(q_w2b, om_b)
    # Newton for moments
    Xd[sv_slice_rvel] = np.dot(P.invJ, M_b - np.cross(om_b, np.dot(P.J, om_b)))
    return Xd
    

def get_forces_and_moments_body(X, U, P):
    # rotors Thrust in body frame
    Fb = [0, 0, -np.sum(U)]
    # Drag
    Dw = -P.Cd*X[sv_slice_vel]
    R_w2b =  pal.rmat_of_quat(X[sv_slice_quat])
    Db = np.dot(R_w2b, Dw)
    
    # Moments of external forces
    Mb = np.sum([np.cross(_p, [0, 0, -_f]) for _p, _f in zip(P.rotor_pos, U)], axis=0)
    # Torques
    Mb[2] += np.sum(P.k*(P.rotor_dir*U))
    return Fb+Db, Mb
    
    
def cont_dyn(X, t, U, P):
    '''
    Continuous-time State Space Representation: Xdot = f_param(t, X, U)
    '''
    Fb, Mb = get_forces_and_moments_body(X, U, P)
    Xd = solid_cont_dyn(X, Fb, Mb, P)
    return Xd


def disc_dyn(Xk, tk, Uk, dt, P):
    '''
    Discrete-time State Space Representation: Xk+1 = f_param(Xk, Uk)
    '''
    _unused, Xkp1 = scipy.integrate.odeint(cont_dyn, Xk, [tk, tk+dt], args=(Uk, P))
    #Xkp1 = Xk + cont_dyn(Xk, tk, Uk, P)*dt # first order
    # normalize quaternion ?
    return Xkp1





def num_jacobian(X, U, P):
  epsilonX = (0.01*np.ones(sv_size)).tolist()
  dX = np.diag(epsilonX)
  A = np.zeros((sv_size, sv_size))
  for i in range(0, sv_size):
    dx = dX[i,:]
    delta_f = cont_dyn(X+dx/2, 0, U, P) - cont_dyn(X-dx/2, 0, U, P)
    delta_f = delta_f / dx[i]
    A[:,i] = delta_f

  epsilonU = (0.1*np.ones(iv_size)).tolist()
  dU = np.diag(epsilonU)
  B = np.zeros((sv_size,iv_size))
  for i in range(0, iv_size):
    du = dU[i,:]
    delta_f = cont_dyn(X, 0, U+du/2, P) - cont_dyn(X, 0, U-du/2, P)
    delta_f = delta_f / du[i]
    B[:,i] = delta_f

  return A,B



def plot(time, X, U=None, figure=None, window_title="Trajectory"):
  figure = ppu.prepare_fig(figure, window_title, (20.48, 10.24))
 
  eulers = np.array([pal.euler_of_quat(_q) for _q in X[:,sv_slice_quat]])
  phi, theta, psi = eulers[:,0], eulers[:,1], eulers[:,2]
  plots = [("$x$",       "m",     0.5, X[:,sv_x]),
           ("$y$",       "m",     0.5, X[:,sv_y]),
           ("$z$",       "m",     0.5, X[:,sv_z]),
           ("$\dot{x}$", "m/s",   0.5, X[:,sv_xd]),
           ("$\dot{y}$", "m/s",   0.5, X[:,sv_yd]),
           ("$\dot{z}$", "m/s",   0.5, X[:,sv_zd]),
           ("$\phi$",    "deg",   0.5, np.rad2deg(phi)),
           ("$\\theta$", "deg",   0.5, np.rad2deg(theta)),
           ("$\psi$",    "deg",   0.5, np.rad2deg(psi)),
           ("$p$",       "deg/s", 0.5, np.rad2deg(X[:,sv_p])),
           ("$q$",       "deg/s", 0.5, np.rad2deg(X[:,sv_q])),
           ("$r$",       "deg/s", 0.5, np.rad2deg(X[:,sv_r])),
  ]
  if U is not None:
    foo = np.empty((len(time))); foo.fill(np.nan)
    plots += [("$U$", "N", 0.1, foo)]
  figure = ppu.plot_in_grid(time, plots, 3, figure, window_title)
  if U is not None:
    ax = plt.subplot(5, 3, 13)
    for i,txt in enumerate(('fr', 'br', 'bl', 'fl')):
      plt.plot(time, U[:,i], label=txt)
    plt.legend()
  return figure
  
  pass
