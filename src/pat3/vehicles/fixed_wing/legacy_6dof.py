#-*- coding: utf-8 -*-
#
# Copyright 2013-2014 Antoine Drouin (poinix@gmail.com)
#
# This file is part of PAT.
#
#    PAT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PAT is distributed in the hpuope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PAT.  If not, see <http://www.gnu.org/licenses/>.
#
import pdb
"""
 This is a 6dof model for a fixed wing vehicle

 legacy code hacked from pat1. \strike{Here to serve as reference}
 I'm hacking that quite a bit :(

"""

import os
import math, numpy as np
import scipy.integrate, scipy.optimize
import matplotlib.pyplot as plt
import importlib

import pat3.utils as p3_u
import pat3.frames as p3_fr
import pat3.plot_utils as p3_pu
import pat3.atmosphere as p3_atm
import pat3.algebra as p3_alg

def load(type, variant, param=None):
    mod_str = "pat.vehicles.{:s}.dynamic_model_{:s}".format(type, variant)
    dm = importlib.import_module(mod_str)
    return dm.DynamicModel(param)


"""
Components of the state vector
"""
sv_x     = 0   # position x axis
sv_y     = 1   # position y axis
sv_z     = 2   # heigh above ground
sv_v     = 3   # airspeed
sv_alpha = 4   # alpha
sv_beta  = 5   # beta
sv_phi   = 6   # roll  (euler, ltp to body)
sv_theta = 7   # pitch (euler, ltp to body)
sv_psi   = 8   # yaw   (euler, ltp to body)
sv_p     = 9   # rotational vel body x
sv_q     = 10  # rotational vel body y
sv_r     = 11  # rotational vel body z
sv_size  = 12

sv_slice_pos   = slice(sv_x,   sv_z+1)
sv_slice_vaero = slice(sv_v,   sv_beta+1)
sv_slice_eul   = slice(sv_phi, sv_psi+1)
sv_slice_rvel  = slice(sv_p,   sv_r+1)


"""
Components of the input vector
"""
iv_da   = 0
iv_de   = 1
iv_dr   = 2
iv_df   = 3
iv_size = 4


# REMOVE this... now in frames
def get_aero_to_body(X):
    """
    computes the aero to body rotation matix
    """
    ca = math.cos(X[sv_alpha]); sa = math.sin(X[sv_alpha]) 
    cb = math.cos(X[sv_beta]);  sb = math.sin(X[sv_beta])
    return np.array([[ca*cb, -ca*sb, -sa],
                     [sb   ,  cb   ,  0.],
                     [sa*cb, -sa*sb,  ca]])

def thrust_of_throttle(throttle, fmax, rho, rhor, nrho, v, vr, nv):
    return throttle*fmax*math.pow((rho/rhor),nrho)*math.pow((v/vr),nv)

def get_f_eng_body2(h, va, U, P):
    """
    return propulsion forces expressed in body frame
    """
    rho = p3_atm.get_rho(h)
    f_engines_body = np.zeros((P.eng_nb, 3))
    for i in range(0, P.eng_nb): # FIXME use list constructor
        thrust = thrust_of_throttle(U[i], P.fmaxs[i], rho, P.rhois[i], P.nrhos[i], va, P.Vis[i], P.nVs[i])
        f_engines_body[i] = np.dot(P.eng_to_body[i], np.array([thrust, 0., 0.]))
    return f_engines_body 

def get_f_eng_body(X, U, P):
    """
    return propulsion forces expressed in body frame
    """
    rho = p3_atm.get_rho(-X[sv_z])
    f_engines_body = np.zeros((P.eng_nb, 3))
    for i in range(0, P.eng_nb):
        thrust = thrust_of_throttle(U[i], P.fmaxs[i], rho, P.rhois[i], P.nrhos[i], X[sv_v], P.Vis[i], P.nVs[i])
        f_engines_body[i] = np.dot(P.eng_to_body[i], np.array([thrust, 0., 0.]))
    return f_engines_body 

def get_f_aero_coef(alpha, beta, rvel, Usfc, P):
    d_alpha = alpha - P.alpha0
    nrvel = rvel*np.array([P.Bref, P.Cref, P.Bref])/2/P.Vref # FIXME va??
    CL = P.CL0 + P.CL_alpha*d_alpha + P.CL_beta*beta + np.dot(P.CL_omega,nrvel) + np.dot(P.CL_sfc,Usfc)
    CD = P.CD0 + P.CD_k1*CL + P.CD_k2*(CL**2) + np.dot(P.CD_sfc,Usfc)
    CY = P.CY_alpha*d_alpha + P.CY_beta*beta + np.dot(P.CY_omega,nrvel) + np.dot(P.CY_sfc,Usfc)
    return [CL, CY, CD]


def get_f_aero_body2(va, alpha, beta, rvel, Usfc, P, Pdyn):
    CL, CY, CD = get_f_aero_coef(alpha, beta, rvel, Usfc, P)
    F_aero_body = Pdyn*P.Sref*np.dot(p3_fr.R_aero_to_body(alpha, beta), [-CD, CY, -CL])
    return F_aero_body


def get_f_aero_body(X, Usfc, P, Pdyn):
    """
    return aerodynamic forces expressed in body frame
    """
    if 0:
        d_alpha = X[sv_alpha] - P.alpha0
        rvel =  X[sv_slice_rvel]*np.array([P.Bref, P.Cref, P.Bref])/2/P.Vref

        CL = P.CL0 + P.CL_alpha*d_alpha + P.CL_beta*X[sv_beta] +\
             np.dot(P.CL_omega,rvel) + np.dot(P.CL_sfc,Usfc)

        CD = P.CD0 + P.CD_k1*CL + P.CD_k2*(CL**2) + np.dot(P.CD_sfc,Usfc)

        CY = P.CY_alpha*d_alpha + P.CY_beta*X[sv_beta] +\
             np.dot(P.CY_omega,rvel) + np.dot(P.CY_sfc,Usfc)
        alpha, beta = X[sv_alpha], X[sv_beta]
        F_aero_body = Pdyn*P.Sref*np.dot(p3_fr.R_aero_to_body(alpha, beta),[-CD, CY, -CL])
    if 1:
        alpha, beta, rvel = X[sv_alpha], X[sv_beta], X[sv_slice_rvel]
        CL, CY, CD = get_f_aero_coef(alpha, beta, rvel, Usfc, P)
        F_aero_body = Pdyn*P.Sref*np.dot(p3_fr.R_aero_to_body(alpha, beta),[-CD, CY, -CL])
    return F_aero_body
    
def get_m_eng_body(f_eng_body, P):
    """
    return propulsion moments expressed in body frame
    """
    m = np.zeros(3)
    for i in range(0, P.eng_nb):
        m += np.cross(P.eng_pos[i], f_eng_body[i])
    return m


def get_m_aero_coef(alpha, beta, rvel, Usfc, P):
    d_alpha = alpha - P.alpha0
    Cl =         P.Cl_alpha*d_alpha + P.Cl_beta*beta +\
         np.dot(P.Cl_omega,rvel) + np.dot(P.Cl_sfc,Usfc)
    Cm = P.Cm0 + P.Cm_alpha*d_alpha + P.Cm_beta*beta +\
         np.dot(P.Cm_omega,rvel) + np.dot(P.Cm_sfc,Usfc)
    Cn =         P.Cn_alpha*d_alpha + P.Cn_beta*beta +\
         np.dot(P.Cn_omega,rvel) + np.dot(P.Cn_sfc,Usfc)
    return Cl, Cm, Cn

def get_m_aero_body2(va, alpha, beta, rvel, Usfc, P, Pdyn):
    Cl, Cm, Cn = get_m_aero_coef(alpha, beta, rvel, Usfc, P)
    return Pdyn*P.Sref*np.array([Cl*P.Bref, Cm*P.Cref, Cn*P.Bref])

def get_m_aero_body(X, Usfc, P, Pdyn):
    """
    return aerodynamic moments expressed in body frame
    """
    d_alpha = X[sv_alpha] - P.alpha0
    rvel =  X[sv_p:sv_r+1]*np.array([P.Bref, P.Cref, P.Bref])/2/P.Vref

    Cl =         P.Cl_alpha*d_alpha + P.Cl_beta*X[sv_beta] +\
         np.dot(P.Cl_omega,rvel) + np.dot(P.Cl_sfc,Usfc)
    Cm = P.Cm0 + P.Cm_alpha*d_alpha + P.Cm_beta*X[sv_beta] +\
         np.dot(P.Cm_omega,rvel) + np.dot(P.Cm_sfc,Usfc)
    Cn =         P.Cn_alpha*d_alpha + P.Cn_beta*X[sv_beta] +\
         np.dot(P.Cn_omega,rvel) + np.dot(P.Cn_sfc,Usfc)

    return Pdyn*P.Sref*np.array([Cl*P.Bref, Cm*P.Cref, Cn*P.Bref])

def dyn(X, t, U, P, atm=None):
    """
    Dynamic model
    """
    rho = p3_atm.get_rho(-X[sv_z])
    Pdyn = 0.5*rho*X[sv_v]**2

    Ueng = U[0:P.eng_nb]                    # engines part of input vector
    Usfc = U[P.eng_nb:P.eng_nb+P.sfc_nb]    # control surfaces part of input vector
    X_pos = X[sv_slice_pos]                 #
    X_rvel_body = X[sv_slice_rvel]          # body rotational velocities
    X_euler = X[sv_slice_eul]               # euler angles                          

    aero_to_body = get_aero_to_body(X)
    earth_to_body = p3_alg.rmat_of_euler(X_euler)
    body_to_earth = earth_to_body.T 
    avel_aero = [X[sv_v], 0., 0.]
    avel_body = np.dot(aero_to_body, avel_aero)
    wvel_earth = (atm.get_wind(X_pos, t) if atm is not None else [0, 0, 0])
    wvel_body = np.dot(earth_to_body, wvel_earth)
    waccel_body = [0, 0, 0]  # np.dot(earth_to_body, waccel_earth)
    ivel_body = avel_body + wvel_body
    
    # Newton for forces in body frame
    f_aero_body = get_f_aero_body(X, Usfc, P, Pdyn)
    f_eng_body  = get_f_eng_body(X, Ueng, P)
    f_weight_body = np.dot(earth_to_body, [0., 0., P.m*P.g])
    forces_body = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body
    
    iaccel_body = 1./P.m*forces_body - np.cross(X_rvel_body, ivel_body)

    # Newton for moments in body frame
    m_aero_body = get_m_aero_body(X, Usfc, P, Pdyn)
    m_eng_body = get_m_eng_body(f_eng_body, P)
    raccel_body = np.dot(P.invI, m_aero_body + m_eng_body - np.cross(X_rvel_body, np.dot(P.I, X_rvel_body)))
    #raccel_body = np.dot(P.invJ, m_aero_body + m_eng_body - np.cross(X_rvel_body, np.dot(P.J, X_rvel_body)))

    Xdot = np.zeros(p3_fr.SixDOFAeroEuler.sv_size)
    Xdot[sv_x:sv_z+1] = np.dot(body_to_earth, ivel_body)

    # CHECK...
    # Air velocity kinematics
    aaccel_body = iaccel_body - waccel_body

    Xdot[sv_v] = np.inner(avel_body, aaccel_body)/X[sv_v]
    #u, v, w = avel_body
    #ud, vd, wd = aaccel_body
    #Xdot[sv_alpha] = (u*wd - w*ud)/(u**2+w**2)
    #Xdot[sv_beta] = (X[sv_v]*vd - v*Xdot[sv_v]) / X[sv_v] / math.sqrt(u**2+w**2)
    (avx, avy, avz), (aax, aay, aaz) = avel_body, aaccel_body
    Xdot[sv_alpha] = (avx*aaz - avz*aax)/(avx**2+avz**2)
    Xdot[sv_beta] = (X[sv_v]*aay - avy*Xdot[sv_v]) / X[sv_v] / math.sqrt(avx**2+aaz**2)
    

    # Euler angles kinematics
    Xdot[sv_phi:sv_psi+1] = p3_alg.euler_derivatives(X_euler, X_rvel_body)
    
    Xdot[sv_p:sv_r+1] = raccel_body

    return Xdot


def trim(P, args, report=True, debug=False):
    """
    Find throttle, elevator  and angle of attack corresponding
    to the given airspeed and and flight path
    """
    flaps = args.get('flaps', 0.)
    if 'gamma' in args:
        va, gamma_e, h = args['va'], args['gamma'], args['h']
        thr_e, ele_e, alpha_e = trim_cst_path(P, h, va, gamma_e, report, debug)
    else:
        va, h, thr_e = args['va'], args['h'], args['throttle']
        gamma_e, ele_e, alpha_e = trim_cst_throttle(P, h, va, thr_e, flaps, report, debug)
        
    Ue = np.zeros(P.input_nb)
    Ue[0:P.eng_nb] = thr_e; Ue[P.eng_nb+iv_de] = ele_e; Ue[P.eng_nb+iv_df] = flaps
    Xe = [0., 0., -h, va, alpha_e, 0., 0., gamma_e+alpha_e, 0., 0., 0., 0.]  
    return Xe, Ue
        
def trim_cst_path(P, h=0., va=12., gamma=0., report=True, debug=False):
    if report:
        print("searching for constant path trajectory with")
        print("  h      {:.1f} m".format(h))
        print("  va     {:.2f} m/s".format(va))
        print("  gamma  {:.2f} deg".format(np.rad2deg(gamma)))

    def err_func(args):
        throttle, elevator, alpha = args
        X=[0., 0., -h, va, alpha, 0., 0.,  gamma+alpha, 0., 0., 0., 0.]
        U = np.zeros(P.input_nb)
        U[0:P.eng_nb] = throttle; U[P.eng_nb+iv_de] = elevator
        Xdot = dyn(X, 0., U, P)
        Xdot_ref = [va*math.cos(gamma), 0., -va*math.sin(gamma), 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        return np.linalg.norm(Xdot - Xdot_ref)

    p0 = [0.2, np.deg2rad(2.), np.deg2rad(0.)]
    thr_e, ele_e, alpha_e = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-9)

    if report:
        print("""result:
  throttle        : {:.2f} %
  elevator        : {:.2f} deg
  angle of attack : {:.2f} deg""".format(100.*thr_e, np.rad2deg(ele_e), np.rad2deg(alpha_e)))

    return thr_e, ele_e, alpha_e


def trim_cst_throttle(P, h=0., va=12., throttle=0., flaps=0., report=True, debug=False):
    if report:
        print("searching for constant throttle trajectory with")
        print("  h      {:.1f} m".format(h))
        print("  va     {:.2f} m/s".format(va))
        print("  throttle  {:.1f} %".format(100*throttle))
        print("  flap  {:.1f} deg".format(np.rad2deg(flaps)))
    def err_func(args):
        gamma, elevator, alpha = args
        X=[0., 0., -h, va, alpha, 0., 0.,  gamma+alpha, 0., 0., 0., 0.]
        U = np.zeros(P.input_nb)
        U[0:P.eng_nb] = throttle; U[P.eng_nb+iv_de] = elevator; U[P.eng_nb+iv_df] = flaps
        Xdot = dyn(X, 0., U, P)
        Xdot_ref = [va*np.cos(gamma), 0., -va*np.sin(gamma), 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        return np.linalg.norm(Xdot - Xdot_ref)

    p0 = [0., np.deg2rad(-2.), np.deg2rad(1.)]
    gamma_e, ele_e, alpha_e = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-9)
    if report:
        print("""result:
  gamma           : {:.2f} deg
  elevator        : {:.2f} deg
  angle of attack : {:.2f} deg""".format(np.rad2deg(gamma_e), np.rad2deg(ele_e), np.rad2deg(alpha_e)))
    return gamma_e, ele_e, alpha_e
        

def trim_banked(P, h=0., va=12., throttle=0., phi=0., report=True, debug=False):
    if report:
        print("searching for banked constant throttle trajectory with")
        print("  roll  {:.1f} deg".format(np.rad2deg(phi)))
        print("  h      {:.1f} m".format(h))
        print("  va     {:.2f} m/s".format(va))
        print("  throttle  {:.1f} %".format(100*throttle))
        R = va**2/np.tan(phi)
        psi = 0. # does not depend on yaw... or?
        phid_ref, thetad_ref, psid_ref = 0., 0., va/R
        #euler_ref, eulerd_ref = [phi, theta, psi], [0, 0, psi_dot]
        #pref, qref, rref = 0, 0, va/R
        print(" psid {} deg/s {}s/rev".format(np.rad2deg(psid_ref), 360./np.rad2deg(psid_ref)))

        
    def err_func(args):
        gamma, elevator, alpha = args
        theta = gamma + alpha
        psi_dot = va/R # not really....
        euler, eulerd = [phi, theta, psi], [0, 0, psi_dot]
        #p, q, r = 0, 0, va/R # FIXME...
        #p, q, r = p3_alg.rvel_of_eulerd(euler, eulerd)
        p, q, r = 0, np.sin(phi)*np.cos(theta)*va/R, np.cos(phi)*np.cos(theta)*va/R
        X=[0., 0., -h, va, alpha, 0., phi,  theta, psi, p, q, r]
        U = np.zeros(P.input_nb)
        U[0:P.eng_nb] = throttle; U[P.eng_nb+iv_de] = elevator
        Xdot = dyn(X, 0., U, P)
        Xdot_ref = [va*math.cos(gamma), 0., -va*math.sin(gamma), 0., 0., 0., phid_ref, thetad_ref, psid_ref, 0, 0, 0]
        pdb.set_trace()
        return np.linalg.norm(Xdot - Xdot_ref)

    p0 = [0., np.deg2rad(-2.), np.deg2rad(1.)]
    gamma_e, ele_e, alpha_e = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-9)
    if report:
        print("""result:
  gamma           : {:.2f} deg
  elevator        : {:.2f} deg
  angle of attack : {:.2f} deg""".format(np.rad2deg(gamma_e), np.rad2deg(ele_e), np.rad2deg(alpha_e)))
    return gamma_e, ele_e, alpha_e
    


import pat3.vehicles.fixed_wing.legacy_parameter
class ParamOld(pat3.vehicles.fixed_wing.legacy_parameter.Param):
    pass
import pat3.vehicles.fixed_wing.simple_6dof_fdm_param
class Param(pat3.vehicles.fixed_wing.simple_6dof_fdm_param.Param):
    pass



#class DynamicModel(BaseDynamicModel):
class DynamicModel:

    sv_x     = sv_x     # position x axis
    sv_y     = sv_y     # position y axis
    sv_z     = sv_z     # heigh above ground
    sv_v     = sv_v     # airspeed
    sv_alpha = sv_alpha # alpha
    sv_beta  = sv_beta  # beta
    sv_phi   = sv_phi   # roll  (euler, ltp to body)
    sv_theta = sv_theta # pitch (euler, ltp to body)
    sv_psi   = sv_psi   # yaw   (euler, ltp to body)
    sv_p     = sv_p     # rotational vel body x
    sv_q     = sv_q     # rotational vel body y
    sv_r     = sv_r     # rotational vel body z
    sv_size  = sv_size

    
    iv_th   = 0    # throttle
    iv_da   = 1    # aileron
    iv_de   = 2    # elevator
    iv_dr   = 3    # rudder
    iv_size = 4
    # hack for multiple engines
    _iv_da = 0
    _iv_de = 1
    _iv_dr = 2
    _iv_df = 3

    # looks broken... dyn = lambda self, X, t, U, P: dyn(X, t, U, self.P)
    trim = lambda self, args=None, report=False, debug=False: trim(self.P, args, report, debug)
    def dyn(self, X, t, U, atm=None):
        return dyn(X, t, U, self.P, atm)

    
    def __init__(self, params=None):
        print("Info: Dynamic fixed wing legacy")
        #BaseDynamicModel.__init__(self)
        if params == None: params = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
        self.P = Param(params)
        #self.P = ParamOld(params)
        self.X = np.zeros(self.sv_size())
        self.X_act = np.zeros(self.input_nb())
        self.t, self.dt = 0., 0.01
        self.T_w2b = np.eye(4)
        self.reset()

    def sv_size(self): return p3_fr.SixDOFAeroEuler.sv_size
 
    def name(self):
        return "Fixed Wing Python Basic ({:s})".format(self.P.name)

    def reset(self, X0=None, t0=0, X_act0=None):
        if X0 is not None: self.X = np.asarray(X0)
        else: self.X = np.array([0., 0., 0., 68., 0., 0., 0., 0., 0., 0., 0., 0.])
        if X_act0 is not None:
            self.X_act = np.asarray(X_act0)
        self.t = t0
        self._update_byproducts()
        return self.X

    def _clip_input(self, U):
        Uclip = np.copy(U)
        # clip input
        max_ail, max_ele, max_rud = np.deg2rad(30), np.deg2rad(25), np.deg2rad(30) 
        Uclip[self.iv_dth()] = np.clip(U[self.iv_dth()], 0., 1.)
        Uclip[self.iv_da()]  = np.clip(U[self.iv_da()], -max_ail, max_ail)
        Uclip[self.iv_de()]  = np.clip(U[self.iv_de()], -max_ele, max_ele)
        Uclip[self.iv_dr()]  = np.clip(U[self.iv_dr()], -max_rud, max_rud)
        return Uclip
        
    def run(self, dt, tf, U, atm):
        if 0: # no actuator dynamics
            self.X_act = self._clip_input(U)
        else:
            act_lambda = -1./np.array([0.1, 0.02, 0.02, 0.02, 0.02])
            self.X_act += dt*((self.X_act-self._clip_input(U))*act_lambda)
       
        foo, self.X = scipy.integrate.odeint(dyn, self.X, [self.t, self.t+dt], args=(self.X_act, self.P, atm))#, hmax=0.001)
        self.t += dt
        self._update_byproducts()
        return self.X

    def _update_byproducts(self):
        #self.T_w2b[:3,:3] = p3_alg.rmat_of_quat(self.X[sv_slice_quat]).T # that is freaking weird....
        self.T_w2b[:3,:3] = p3_alg.rmat_of_euler(self.X[sv_slice_eul]).T
        self.T_w2b[:3,3] = self.X[sv_slice_pos]
        
    def param(self):
        return str(self.P)

    def iv_dth(self):
        if self.P.eng_nb>1: return range(0,self.P.eng_nb)
        else: return 0
    def iv_da(self): return self.P.eng_nb + DynamicModel._iv_da
    def iv_de(self): return self.P.eng_nb + DynamicModel._iv_de
    def iv_dr(self): return self.P.eng_nb + DynamicModel._iv_dr
    def iv_df(self): return self.P.eng_nb + DynamicModel._iv_df
    def input_nb(self): return self.P.input_nb

    def thrust_of_throttle(self, throttle, i, h, v):
        rho = p3_atm.get_rho(h)
        return thrust_of_throttle(throttle, self.P.fmaxs[i], 
                                  rho, self.P.rhois[i], self.P.nrhos[i], 
                                  v, self.P.Vis[i], self.P.nVs[i])

    def state_six_dof_euclidian_euler(self, X, atm=None):
        return p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(X, atm)
    
    def get_jacobian(self, Xe, Ue):
        A,B = p3_u.num_jacobian(Xe, Ue, self.P, dyn)
        return A, B 

    def state_str(self, X=None):
        return p3_fr.SixDOFAeroEuler.state_str(X if X is not None else self.X)

    def plot_trajectory(self, time, X, U=None, figure=None, window_title="Trajectory", legend=None, filename=None): 
        plot_trajectory(time, X, U, figure, window_title, legend, filename)

    def plot_trajectory_ee(self, time, Xae, U=None, figure=None, window_title="Trajectory", legend=None, filename=None, atm=None): 
        Xee = np.array([p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(_X, atm) for _X in Xae])
        plot_trajectory_ee(time, Xee, U, figure, window_title, legend, filename)



class DynamicModel_ee:

    ss = p3_fr.SixDOFEuclidianEuler
    
    def __init__(self, params_filesname=None):
        print("Info: Dynamic fixed wing legacy ee")
        if params_filesname == None: params_filesname = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
        self.P = Param(params_filesname)
        self.T_w2b = np.eye(4)

    def dyn(self, X, t, U, atm):

        Ueng = U[self.P.u_slice_eng()]                       # engines part of input vector
        Usfc = U[self.P.u_slice_sfc()]                       # control surfaces part of input vector
        X_pos = X[self.ss.sv_slice_pos]                      #
        X_rvel_body = X[self.ss.sv_slice_rvel]               # body rotational velocities
        X_euler = X[self.ss.sv_slice_eul]                    # euler angles                          

        earth_to_body_R = p3_alg.rmat_of_euler(X_euler)
        va, alpha, beta = p3_fr.vel_world_to_aero(X_pos, X[self.ss.sv_slice_vel], X_euler, atm)
        rho = p3_atm.get_rho(-X[self.ss.sv_z])
        Pdyn = 0.5*rho*va**2
        
        Xdot = np.zeros(self.ss.sv_size)
        # Position kinematics
        Xdot[self.ss.sv_slice_pos] = X[self.ss.sv_slice_vel]
        # Newton for forces in world frame
        f_aero_body = get_f_aero_body2(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
        f_eng_body  =  get_f_eng_body2(-X[self.ss.sv_z], va, Ueng, self.P)
        f_weight_body = np.dot(earth_to_body_R, [0., 0., self.P.m*self.P.g])
        forces_body = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body
        forces_world = np.dot(earth_to_body_R.T, forces_body)
        Xdot[self.ss.sv_slice_vel] = 1./self.P.m * forces_world
        # Orientation kinematics
        Xdot[self.ss.sv_slice_eul] = p3_alg.euler_derivatives(X_euler, X_rvel_body)
        # Newton for moments in body frame
        m_aero_body = get_m_aero_body2(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
        m_eng_body = get_m_eng_body(f_eng_body, self.P)
        raccel_body = np.dot(self.P.invI, m_aero_body + m_eng_body - np.cross(X_rvel_body, np.dot(self.P.I, X_rvel_body)))
        Xdot[self.ss.sv_slice_rvel] = raccel_body

        return Xdot

    def reset(self, X0=None, t0=0, X_act0=None):
        if X0 is not None: self.X = np.asarray(X0)
        else: self.X = np.array([0., 0., 0., 68., 0., 0., 0., 0., 0., 0., 0., 0.])
        if X_act0 is not None:
            self.X_act = np.asarray(X_act0)
        self.t = t0
        self._update_byproducts()
        return self.X

    def run(self, dt, tf, U, atm):
        if 0: # no actuator dynamics
            self.X_act = self._clip_input(U)
        else:
            act_lambda = -1./np.array([0.1, 0.02, 0.02, 0.02, 0.02])
            self.X_act += dt*((self.X_act-self._clip_input(U))*act_lambda)
       
        foo, self.X = scipy.integrate.odeint(dyn, self.X, [self.t, self.t+dt], args=(self.X_act, self.P, atm))#, hmax=0.001)
        self.t += dt
        self._update_byproducts()
        return self.X

    def _clip_input(self, U):
        Uclip = np.copy(U)
        # clip input
        max_ail, max_ele, max_rud = np.deg2rad(30), np.deg2rad(25), np.deg2rad(30) 
        Uclip[self.iv_dth()] = np.clip(U[self.iv_dth()], 0., 1.)
        Uclip[self.iv_da()]  = np.clip(U[self.iv_da()], -max_ail, max_ail)
        Uclip[self.iv_de()]  = np.clip(U[self.iv_de()], -max_ele, max_ele)
        Uclip[self.iv_dr()]  = np.clip(U[self.iv_dr()], -max_rud, max_rud)
        return Uclip

    def _update_byproducts(self):
        self.T_w2b[:3,:3] = p3_alg.rmat_of_euler(self.X[self.ss.sv_slice_eul]).T
        self.T_w2b[:3,3] = self.X[self.ss.sv_slice_pos]
 
    def trim(self, args=None, report=False, debug=False):

        if 'gamma' in args:
            va, gamma_e, h = args['va'], args['gamma'], args['h']
            thr_e, ele_e, alpha_e = self.trim_cst_path(h, va, gamma_e, report, debug)
        
        Ue = np.zeros(self.P.u_size())
        Ue[self.P.u_slice_eng()], Ue[self.P.u_elevator()] = thr_e, ele_e
        xd, yd, zd = va*math.cos(gamma_e), 0., -va*math.sin(gamma_e)
        Xe = [0., 0., -h, xd, yd, zd, 0., gamma_e+alpha_e, 0., 0., 0., 0.]  
        return Xe, Ue


    def trim_cst_path(self, h, va, gamma, report, debug, atm=None):
        if report:
            print("searching for constant path trajectory with")
            print("  h      {:.1f} m".format(h))
            print("  va     {:.2f} m/s".format(va))
            print("  gamma  {:.2f} deg".format(np.rad2deg(gamma)))

        wvel_world = (atm.get_wind([0, 0, -h], t=0) if atm is not None else [0, 0, 0])
        
        
        def err_func(args):
            throttle, elevator, alpha = args
            xd, yd, zd = va*math.cos(gamma), 0., -va*math.sin(gamma)
            phi, theta, psi = 0.,  gamma+alpha, 0.
            X=[0., 0., -h, xd, yd, zd, phi, theta, psi, 0., 0., 0.]
            U = np.zeros(self.P.input_nb)
            U[self.P.u_slice_eng()], U[self.P.u_elevator()] = throttle, elevator
            Xdot = self.dyn(X, 0., U, atm)
            Xdot_ref = [va*math.cos(gamma), 0., -va*math.sin(gamma), 0., 0., 0., 0., 0., 0., 0., 0., 0.]
            return np.linalg.norm(Xdot - Xdot_ref)

        p0 = [0.2, np.deg2rad(2.), np.deg2rad(0.)]
        thr_e, ele_e, alpha_e = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-9)

        if report:
            print("""result:
    throttle        : {:.2f} %
    elevator        : {:.2f} deg
    angle of attack : {:.2f} deg""".format(100.*thr_e, np.rad2deg(ele_e), np.rad2deg(alpha_e)))

        return thr_e, ele_e, alpha_e

    def input_nb(self): return self.P.input_nb
    def sv_size(self): return self.ss.sv_size
    def iv_dth(self):
        if self.P.eng_nb>1: return range(0,self.P.eng_nb)
        else: return 0
    def iv_da(self): return self.P.eng_nb + DynamicModel._iv_da
    def iv_de(self): return self.P.eng_nb + DynamicModel._iv_de
    def iv_dr(self): return self.P.eng_nb + DynamicModel._iv_dr
    def iv_df(self): return self.P.eng_nb + DynamicModel._iv_df
 
    def plot_trajectory_ee(self, time, Xee, U=None, figure=None, window_title="Trajectory", legend=None, filename=None): 
        plot_trajectory_ee(time, Xee, U, figure, window_title, legend, filename)

#
# Some plotting functions
#
def plot_trajectory(time, X, U=None, figure=None, window_title="Trajectory",
                    legend=None, filename=None):

        margins=(0.04, 0.05, 0.98, 0.96, 0.20, 0.34)
        figure = p3_pu.prepare_fig(figure, window_title, figsize=(20.48, 10.24), margins=margins)

        plots = [("x", "m", X[:,sv_x]), ("y", "m", X[:,sv_y]), ("z", "m", X[:,sv_z]),
                 ("v",         "m/s",               X[:,sv_v]),
                 ("$\\alpha$", "deg", np.rad2deg(X[:,sv_alpha])),
                 ("$\\beta$",  "deg", np.rad2deg(X[:,sv_beta])),
                 ("$\phi$",    "deg", np.rad2deg(X[:,sv_phi])),
                 ("$\\theta$", "deg", np.rad2deg(X[:,sv_theta])),
                 ("$\\psi$",   "deg", np.rad2deg(X[:,sv_psi])),
                 ("$p$",     "deg/s", np.rad2deg(X[:,sv_p])),
                 ("$q$",     "deg/s", np.rad2deg(X[:,sv_q])),
                 ("$r$",     "deg/s", np.rad2deg(X[:,sv_r]))]

        nrow = 5 if U is not None else 4
        for i, (title, ylab, data) in enumerate(plots):
            ax = plt.subplot(nrow, 3, i+1)
            plt.plot(time, data)
            p3_pu.decorate(ax, title=title, ylab=ylab)
            
        if legend!=None:
            plt.legend(legend, loc='best')

        for i, min_yspan in enumerate([1., 1., 1.,  1., 1., 1.,  1., 1., 1.,  1., 1., 1.]):
            p3_pu.ensure_yspan(plt.subplot(nrow,3,i+1), min_yspan)

        if U is not None:
            ax = figure.add_subplot(5, 3, 13)
            ax.plot(time, 100*U[:, 0])
            p3_pu.decorate(ax, title="$d_{th}$", ylab="%")
            ax = figure.add_subplot(5, 3, 14)
            ax.plot(time, np.rad2deg(U[:, iv_da+1]))
            p3_pu.decorate(ax, title="$d_a$", ylab="deg", min_yspan=1.)
            ax = figure.add_subplot(5, 3, 15)
            ax.plot(time, np.rad2deg(U[:, iv_de+1]))
            p3_pu.decorate(ax, title="$d_e$", ylab="deg", min_yspan=1.)
    
        return figure
        
# FIXME, move that elsewhere
def plot_trajectory_ee(time, Xee, U=None, figure=None, window_title="Trajectory",
                    legend=None, filename=None):
    margins=(0.04, 0.05, 0.98, 0.96, 0.20, 0.34)
    figure = p3_pu.prepare_fig(figure, window_title, figsize=(20.48, 10.24), margins=margins)
    plots = [("x", "m", Xee[:,sv_x]), ("y", "m", Xee[:,sv_y]), ("z", "m", Xee[:,sv_z]),
             ("$\dot{x}$",  "m/s", Xee[:, p3_fr.SixDOFEuclidianEuler.sv_xd]),
             ("$\dot{y}$",  "m/s", Xee[:, p3_fr.SixDOFEuclidianEuler.sv_yd]),
             ("$\dot{z}$",  "m/s", Xee[:, p3_fr.SixDOFEuclidianEuler.sv_zd]),
             ("$\phi$",     "deg",   np.rad2deg(Xee[:,sv_phi])),
             ("$\\theta$",  "deg",   np.rad2deg(Xee[:,sv_theta])),
             ("$\\psi$",    "deg",   np.rad2deg(Xee[:,sv_psi])),
             ("$p$",        "deg/s", np.rad2deg(Xee[:,sv_p])),
             ("$q$",        "deg/s", np.rad2deg(Xee[:,sv_q])),
             ("$r$",        "deg/s", np.rad2deg(Xee[:,sv_r]))]

    nrow = 5 if U is not None else 4
    for i, (title, ylab, data) in enumerate(plots):
        ax = plt.subplot(nrow, 3, i+1)
        plt.plot(time, data)
        p3_pu.decorate(ax, title=title, ylab=ylab)
            
    if legend!=None:
        plt.legend(legend, loc='best')

    for i, min_yspan in enumerate([1., 1., 1.,  1., 1., 1.,  1., 1., 1.,  1., 1., 1.]):
        p3_pu.ensure_yspan(plt.subplot(nrow,3,i+1), min_yspan)

    if U is not None:
        ax = figure.add_subplot(5, 3, 13)
        ax.plot(time, 100*U[:, 0])
        p3_pu.decorate(ax, title="$d_{th}$", ylab="%", min_yspan=0.01)
        ax = figure.add_subplot(5, 3, 14)
        ax.plot(time, np.rad2deg(U[:, iv_da+1]))
        p3_pu.decorate(ax, title="$d_a$", ylab="deg", min_yspan=1.)
        ax = figure.add_subplot(5, 3, 15)
        ax.plot(time, np.rad2deg(U[:, iv_de+1]))
        p3_pu.decorate(ax, title="$d_e$", ylab="deg", min_yspan=1.)
    
    return figure
