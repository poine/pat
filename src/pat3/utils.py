import sys, os
import numpy as np, math
import threading


import pat3.frames as p3_fr
import pdb

def set_rot(T, R): T[:3,:3]=R
def set_trans(T, t): T[:3,3]=t
def get_rot(T): return T[:3,:3]
def get_trans(T): return T[:3,3]

def norm_mpi_pi(v): return ( v + np.pi) % (2 * np.pi ) - np.pi
def norm_0_2pi(v): return ( v + np.pi) % (2 * np.pi )


# we assume this file is pat_dir/src/pat3/utils.py
def pat_dir():
    dirname, filename = os.path.split(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(dirname, '../..'))
    
def step(t, a=-1., p=10., dt=0.): return a if math.fmod(t+dt, p) > p/2 else -a
def step_vec(t, a=-1., p=10., dt=0.): return np.array([step(_t, a, p, dt) for _t in t])

class Sim:

    def __init__(self, fdm, ctl, atm):
        self.fdm, self.ctl, self.atm = fdm, ctl, atm
        self.Yc = np.zeros(ctl.spv_size)
        self.max_n_X = 5000
        self.Xs  = []  # stores original model state over time
        self.Us  = []  # stores model input over time
        self.Xees = [] # stores model state as euclidian/euler
        
    def reset(self, t0, X0=None, U0=None):
        if X0 is None:
            X0, U0 = self.fdm.trim()
        self.Xs = [X0]
        self.Us = [U0]
        self.Xees = [self.fdm.state_as_six_dof_euclidian_euler(X0)]
        return self.fdm.reset(X0, t0, U0)
        
    def run(self, t1):
        Xee = self.fdm.state_as_six_dof_euclidian_euler(self.fdm.X, self.atm)
        U = self.ctl.get(self.fdm.t, self.fdm.X, Xee, self.Yc) # in case we don't run the fdm
        #pdb.set_trace()
        #print('{} sim.run to {:.3f} (orig {:.3f})'.format(threading.currentThread().getName(), t1, self.fdm.t))
        #if t1 >= self.fdm.t+self.fdm.dt:
        #    print('run sim from {:.3f} to {:.3f}'.format(self.fdm.t, t1))
        #else:
        #    print('not running sim at {:.3f} (next end scheduled to {:.3f}'.format(t1, self.fdm.t+self.fdm.dt))
        while t1 - self.fdm.t >= self.fdm.dt:#0:#self.fdm.dt:
            #print(' compute control at {:.3f}'.format(self.fdm.t))
            Xee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(self.fdm.X, self.atm)

            U = self.ctl.get(self.fdm.t, self.fdm.X, Xee, self.Yc)
            #U =  self.Us[0]
            
            #print(' run fdm from {:.3f} to {:.3f}'.format(self.fdm.t, self.fdm.t+self.fdm.dt))
            self.fdm.run(self.fdm.dt, self.fdm.t+self.fdm.dt, U, self.atm)
            self.Xs.append(self.fdm.X)
            self.Us.append(U)
            self.Xees.append(self.fdm.state_as_six_dof_euclidian_euler(self.fdm.X, self.atm))
            if len(self.Xs) > self.max_n_X:
                del self.Xs[0]
                del self.Us[0]
                del self.Xees[0]
        return U, self.fdm.X
        



'''
  Linear reference models

 Adaptive Control with a Nested Saturation Reference Model
 https://pdfs.semanticscholar.org/8b5c/2718be69a651dc3233a934ac05f5edda7ffd.pdf
'''
class LinRef:
    ''' Nested Saturation Linear Reference Model (with first order integration)'''
    def __init__(self, K, sats=None):
        '''
        K: coefficients of the caracteristic polynomial, in ascending powers order,
              highest order ommited (normalized to -1)
        sats: saturations for each order in ascending order
        '''
        self.K = K; self.order = len(K)
        self.sats = sats
        if self.sats is not None:
            self.M = np.array(sats)
            self.M[0:-1] *= K[1:]
            for i in range(0, len(self.M)-1):
                self.M[len(self.M)-2-i] /= np.prod(self.M[len(self.M)-1-i:])
            self.CM = np.cumprod(self.M[::-1])[::-1]
            #print('M', self.M, 'CM', self.CM)
        self.X = np.zeros(self.order+1)

    def run(self, dt, sp):
        self.X[:self.order] += self.X[1:self.order+1]*dt
        e =  np.array(self.X[:self.order]); e[0] -= sp
        if self.sats is None:
            self.X[self.order] = np.sum(e*self.K)
        else:
            self.X[self.order] = 0
            for i in range(0, self.order):
                self.X[self.order] = self.M[i]*np.clip(self.K[i]/self.CM[i]*e[i] + self.X[self.order], -1., 1.)
        return self.X

    def poles(self):
        return np.roots(np.insert(np.array(self.K[::-1]), 0, -1))

    def reset(self, X0=None):
        if X0 is None: X0 = np.zeros(self.order+1)
        self.X = X0


class FirstOrdLinRef(LinRef):
    def __init__(self, tau):
        LinRef.__init__(self, [-1/tau])

class SecOrdLinRef(LinRef):
    def __init__(self, omega, xi, sats=None):
        LinRef.__init__(self, [-omega**2, -2*xi*omega], sats)

def lamdba_to_omxi(l1, l2): return 
def omxi_to_lambda(om, xi): re, im = -xi*om, np.sqrt(1-xi**2)*om; return np.complex(re,im), np.complex(re, -im)
        
"""
Naive numerical differentiation
"""
def num_jacobian(X, U, P, dyn):
    s_size = len(X)
    i_size = len(U)
    epsilonX = (0.1*np.ones(s_size)).tolist()
    dX = np.diag(epsilonX)
    A = np.zeros((s_size, s_size))
    for i in range(0, s_size):
        dx = dX[i,:]
        delta_f = dyn(X+dx/2, 0, U, P) - dyn(X-dx/2, 0, U, P)
        delta_f = delta_f / dx[i]
        A[:,i] = delta_f

    epsilonU = (0.1*np.ones(i_size)).tolist()
    dU = np.diag(epsilonU)
    B = np.zeros((s_size,i_size))
    for i in range(0, i_size):
        du = dU[i,:]
        delta_f = dyn(X, 0, U+du/2, P) - dyn(X, 0, U-du/2, P)
        delta_f = delta_f / du[i]
        B[:,i] = delta_f

    return A,B
