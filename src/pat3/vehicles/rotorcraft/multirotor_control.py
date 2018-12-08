#! /usr/bin/env python
#-*- coding: utf-8 -*-

import math, numpy as np
import pdb
import pat3.algebra as pal
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm

# /media/mint17/home/poine/dissertation/exemples/guidage_d_un_drone_a_poussee_vectorielle/control.py
# /media/mint17/home/poine/dissertation/obsolete/guidage_quad_multirotors.tex

iv_z, iv_qi, iv_qx, iv_qy, iv_qz = range(5)
iv_size = 5

def step(t, a=-1., p=10., dt=0.): return a if math.fmod(t+dt, p) > p/2 else -a

class StepZInput:
    def get(self, t): return [step(t), 1, 0, 0, 0]

class SinZInput:
    def get(self, t): return [0.5*np.sin(t), 1, 0, 0, 0]

class StepEulerInput:
    def __init__(self, _i, _a=np.deg2rad(1.), p=10, dt=5):
        self._i, self._a, self._p, self.dt = _i, _a, p, dt
        
    def get(self, t):
        eu = [0, 0, 0]
        eu[self._i] = step(t, self._a, self._p, self.dt)
        zc, qc = [0], pal.quat_of_euler(eu)
        return np.append(zc, qc)

class CstInput:
    def __init__(self, z, eu):
        self.z = z
        self.q = pal.quat_of_euler(eu)
        
    def get(self, t):
        return np.append(self.z, self.q)

    

class FooController:

    def __init__(self, fdm):
        self.Xe, self.Ue = fdm.trim()
        self.z_ctl = Zctl(self.Ue, self.Xe)
        self.att_ctl = AttCtl(fdm.P)
        

    def get(self, X, Yc):
        zc, qc = Yc[0], Yc[1:]
        Uz = self.z_ctl.run(X, zc)
        Upqr = self.att_ctl.run(X, qc)
        H = np.array([[0.25, -1,  1, -1],
                      [0.25, -1, -1,  1],
                      [0.25,  1, -1, -1],
                      [0.25,  1,  1,  1]])
        U = np.dot(H, np.hstack((Uz, Upqr)))
        #pdb.set_trace()
        return U


class Zctl:
    def __init__(self, Ue, Xe):
        self.Fez, self.Xe = np.sum(Ue), Xe
        self.K = [[0, 0, -1,  0, 0, -1,  0, 0, 0, 0,  0, 0, 0]]
        self.H = np.array([-1])

    def run(self, X, zc):
        Fbz = self.Fez - np.dot(self.K, X-self.Xe) + np.dot(self.H, [zc])
        return Fbz


class AttCtl:
    def __init__(self, P):
        self.P = P           # dynamic model parameters
        self.ref = AttRef()
        self.omega = np.array([20., 20., 15.])
        self.xi = np.array([0.7, 0.7, 0.7])
            
    def run(self, X, qref=[1, 0, 0, 0], rvel_ref=[0, 0, 0]):
        # error quaternion
        err_quat = pal.quat_inv_comp(X[fdm.sv_slice_quat], qref)
        err_quat = pal.quat_wrap_shortest(err_quat)
        # rotational velocities
        delta_rvel = X[fdm.sv_slice_rvel] - rvel_ref
        # rotational acceleration
        racc = -2*self.omega*self.xi*delta_rvel + self.omega**2 * err_quat[pal.q_x:pal.q_z+1]
        #
        Jxx, Jyy, Jzz = np.diag(self.P.J)
        tmp = np.array([(Jzz-Jyy)/Jxx*X[fdm.sv_q]*X[fdm.sv_r],
                        (Jxx-Jzz)/Jyy*X[fdm.sv_p]*X[fdm.sv_r],
                        (Jyy-Jxx)/Jzz*X[fdm.sv_p]*X[fdm.sv_q]]);
        # inertia
        J = np.array([Jxx/self.P.l, Jyy/self.P.l, Jzz/self.P.k])
        Upqr = J * ( racc  + tmp )
        return Upqr


        
class AttRef:
    def __init__(self):
        self.q = pal.quat_null_ixyz()


    def run(self, pqr_sp, dt):
        self.q = pal.quat_integrate(self.q, pqr_sp, dt)
        self.om = pqr_sp
