#! /usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np, math, matplotlib.pyplot as plt
import scipy.integrate
from io import StringIO
import xml.etree.ElementTree as ET
import logging, ast

"""
 This is a 6dof model for a fixed wing vehicle
"""
from pat3.vehicles.fixed_wing.simple_6dof_fdm_param import Param
import pat3.dynamics as pat_dyn
import pat3.algebra as pal
import pat3.flat_world as pfw
import pat3.atmosphere as patm
import pat3.plot_utils as ppu
LOG = logging.getLogger(__name__)
import pdb

def thrust_of_throttle(throttle, fmax, rho, rhor, nrho, v, vr, nv):
    "  Propulsion Model "
    return throttle*fmax*math.pow((rho/rhor),nrho)*math.pow((v/vr),nv)

def get_f_eng_body(rho, va, Ueng, P):
    " compute propulsion forces expressed in body frame "
    f_engines_body = np.zeros((P.eng_nb, 3))
    for i in range(0, P.eng_nb):
        thrust = thrust_of_throttle(Ueng[i], P.fmaxs[i], rho, P.rhois[i], P.nrhos[i], va, P.Vis[i], P.nVs[i])
        f_engines_body[i] = np.dot(P.eng_to_body[i], np.array([thrust, 0., 0.]))
    return f_engines_body

def get_f_aero_body(alpha, beta, omega, Pdyn, Usfc, P):
    " compute aerodynamic forces expressed in body frame "
    d_alpha = alpha - P.alpha0
    #NEW
    #n_omega =  omega#*np.array([P.Bref, P.Cref, P.Bref]) #/2/P.Vref
    # OLD
    n_omega =  omega*np.array([P.Bref, P.Cref, P.Bref])/2/P.Vref
    #
    CL = P.CL0 + P.CL_alpha*d_alpha + P.CL_beta*beta +\
         np.dot(P.CL_omega, n_omega) + np.dot(P.CL_sfc,Usfc)

    CD = P.CD0 + P.CD_k1*CL + P.CD_k2*(CL**2) + np.dot(P.CD_sfc,Usfc)

    CY = P.CY_alpha*d_alpha + P.CY_beta*beta +\
         np.dot(P.CY_omega, n_omega) + np.dot(P.CY_sfc,Usfc)
    
    return Pdyn*P.Sref*np.dot(pfw.aero_to_body(alpha, beta),[-CD, CY, -CL])


def get_forces_body(rho, va, alpha, beta, eulers, omega, Pdyn, Ueng, Usfc, P):
    " compute sum of external forces expressed in body frame "
    f_aero_body = get_f_aero_body(alpha, beta, omega, Pdyn, Usfc, P)
    f_eng_body  = get_f_eng_body(rho, va, Ueng, P)
    if 0:
        f_weight_body = np.dot(pfw.earth_to_body(eulers), [0., 0., P.m*P.g])
        Fb = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body
    else:
        Fb = f_aero_body + np.sum(f_eng_body, axis=0) 
    return Fb, f_eng_body
    
def get_m_eng_body(f_eng_body, P):
    " propulsion moments expressed in body frame "
    m = np.zeros(3)
    for i in range(0, P.eng_nb):
        m += np.cross(P.eng_pos[i], f_eng_body[i])
    return m

def get_m_aero_body(alpha, beta, euler, omega, Pdyn, Usfc, P):
    " compute aerodynamic moments expressed in body frame "
    d_alpha = alpha - P.alpha0
    n_omega =  omega #_aero*np.array([P.Bref, P.Cref, P.Bref])/2/P.Vref

    Cl =         P.Cl_alpha*d_alpha + P.Cl_beta*beta +\
          np.dot(P.Cl_omega,n_omega) + np.dot(P.Cl_sfc,Usfc)
    Cm = P.Cm0 + P.Cm_alpha*d_alpha + P.Cm_beta*beta +\
          np.dot(P.Cm_omega,n_omega) + np.dot(P.Cm_sfc,Usfc)
    Cn =         P.Cn_alpha*d_alpha + P.Cn_beta*beta +\
          np.dot(P.Cn_omega,n_omega) + np.dot(P.Cn_sfc,Usfc)
    
    m_body = Pdyn*P.Sref*np.array([Cl*P.Bref, Cm*P.Cref, Cn*P.Bref]) 
    #m_body =  np.dot(fr.aero_to_body(alpha, beta), m_aero)

    return  m_body 

def get_moments_body(alpha, beta, euler, omega, Pdyn, f_eng_body, Usfc, P):
    " compute sum of external moments expressed in body frame "
    if 0:
        m_aero_body = get_m_aero_body(alpha, beta, euler, omega, Pdyn, Usfc, P)
        m_eng_body = get_m_eng_body(f_eng_body, P)
        return  m_aero_body + m_eng_body
    else:
        return  get_m_aero_body(alpha, beta, euler, omega, Pdyn, Usfc, P)


 

''' 
  Default actuator allocation used for triming.
  We assume throttle is even on all engines and elevator acts on the second surface.
'''
def Ueng_of_throttle(throttle, P): return [throttle for i in range(0, P.eng_nb)] 
def Usfc_of_elevator(ele, P): return [0 if i!=1 else ele for i in range(0, P.sfc_nb)]


'''
   Flight Dynamic Model for a fixed wing vehicle
'''
class FDM(pat_dyn.SolidFDM):

    def __init__(self, P):
        pat_dyn.SolidFDM.__init__(self, P)
        self.iv_size = P.input_nb
        #print(self.P)
    
    def cont_dyn(self, X, t, U):
        #Fb = np.array([0, 0, 0])
        #Mb = np.array([0, 0, 0])

        Ueng, Usfc = U[:self.P.eng_nb], U[self.P.eng_nb:]   # engines and control surfaces parts of input vector
        #LOG.info(" engine {} sfc {}".format(Ueng, Usfc))
        euler = pal.euler_of_quat(X[self.sv_slice_quat])
        omega = X[self.sv_slice_rvel]
        R_b2w = pal.rmat_of_quat(X[pat_dyn.SolidFDM.sv_slice_quat])
        #pdb.set_trace()
        vi_e = X[self.sv_slice_vel]
        w_e = [0, 0, 0]
        va, alpha, beta = pfw.get_va_alpha_beta(vi_e, w_e, euler) # fixme: make a fonction that takes vi_b
        rho = patm.get_rho(-X[self.sv_z])
        Pdyn = 0.5*rho*va**2
        
        Fb, f_eng_body = get_forces_body(rho, va, alpha, beta, euler, omega, Pdyn, Ueng, Usfc, self.P)
        Mb = get_moments_body(alpha, beta, euler, omega, Pdyn, f_eng_body, Usfc, self.P)
        print(Usfc, Mb)
        # try:
        #     return ps.s1_dyn(X, t, Fb, Mb, P.m, P.I, P.invI, P.inv_mamat)
        # except AttributeError:
        #     return ps.s1_dyn(X, t, Fb, Mb, P.m, P.I, P.invI, None)
        return pat_dyn.SolidFDM.cont_dyn(self, X, t, np.concatenate((Fb, Mb)))
        
   



        
    def trim(self, args={}, UengFun=Ueng_of_throttle, UsfcFun=Usfc_of_elevator, debug=False):
        #Xe = np.zeros(self.sv_size); Xe[self.sv_qi] = 1.; Xe[self.sv_xd] = 15.
        #Ue = np.zeros(self.iv_size)
        
        va, gamma, h = args.get('va', self.P.Vref), args.get('gamma', 0.), args.get('h', 0.)
        if debug:
            print("searching for constant path trajectory with")
            print("  h      {:f} m".format(h))
            print("  va     {:f} m/s".format(va))
            print("  gamma  {:f} deg".format(np.rad2deg(gamma)))
        LOG.info(" Running trim for va {} gamma {} h {}".format(va, gamma, h))
        rho = patm.get_rho(h)
        Pdyn, beta, omega = 0.5*rho*va**2, 0, [0, 0, 0]
        def err_func(args):
            ele, throttle, alpha = args
            Ueng, Usfc = UengFun(throttle, self.P), UsfcFun(ele, self.P)
            eulers = [0, gamma+alpha, 0]
            Fb, f_eng_body = get_forces_body(rho, va, alpha, beta, eulers, omega, Pdyn, Ueng, Usfc, self.P)
            f_weight_body = np.dot(pfw.earth_to_body(eulers), [0., 0., self.P.m*self.P.g])
            Fb = Fb + f_weight_body
            Mb = get_moments_body(alpha, beta, eulers, omega, Pdyn, f_eng_body, Usfc, self.P)
            return [Fb[0], Fb[2], Mb[1]]
        sol = scipy.optimize.root(err_func, [0., 0.5, 0.], method='hybr')
        ele, throttle, alpha = sol.x
        #print sol
        LOG.info(" elevator {:.1f} deg throttle {:.1f} % alpha {:.1f} deg".format(np.rad2deg(ele), 100*throttle, np.rad2deg(alpha)))
        if debug:
            print("""result:
  throttle        : {:f} %
  elevator        : {:f} deg
  angle of attack : {:f} deg""".format(100.*throttle, np.rad2deg(ele), np.rad2deg(alpha)) )

        #Xe = [0, 0, -h, va*math.cos(alpha), 0, va*math.sin(alpha), 0, gamma+alpha, 0, 0, 0, 0]
        q = pal.quat_of_euler([0, gamma+alpha, 0])
        Xe = [0, 0, -h, va*math.cos(gamma), 0, va*math.sin(gamma), q[0], q[1], q[2], q[3], 0, 0, 0]
        Ue = UengFun(throttle, self.P) + UsfcFun(ele, self.P)
        return Xe, Ue

    # def plot(self, time, X, U=None, figure=None, window_title="FW Trajectory"):
    #     Usolid = np.zeros((len(time), pat_dyn.SolidFDM.iv_size))
    #     # Usolid[:,SolidFDM.iv_fzb] = U[:,self.iv_fz]
    #     # Usolid[:,SolidFDM.iv_mxb] = U[:,self.iv_mx]
    #     # Usolid[:,SolidFDM.iv_myb] = U[:,self.iv_my]
    #     # Usolid[:,SolidFDM.iv_mzb] = U[:,self.iv_mz]
    #     return pat_dyn.SolidFDM.plot(self, time, X, Usolid, figure, window_title)




class FDM3(pat_dyn.SolidDM3):

    def __init__(self, P):
        pat_dyn.SolidDM3.__init__(self, P)
        self.iv_size = P.input_nb
        #print(self.P)

    def input_nb(self): return self.P.input_nb
        
    def cont_dyn(self, X, t, U):
        Ueng, Usfc = U[:self.P.eng_nb], U[self.P.eng_nb:]   # engines and control surfaces parts of input vector

        rho = patm.get_rho(-X[self.sv_z])
        Pdyn = 0.5*rho*X[self.sv_v]**2
        (va, alpha, beta), euler, omega = X[self.sv_slice_vaero], X[self.sv_slice_eul], X[self.sv_slice_rvel]
    
        Fb, f_eng_body = get_forces_body(rho, va, alpha, beta, euler, omega, Pdyn, Ueng, Usfc, self.P)
        Mb = get_moments_body(alpha, beta, euler, omega, Pdyn, f_eng_body, Usfc, self.P)
        #pdb.set_trace()
        #print Usfc, Mb
        return pat_dyn.SolidDM3.cont_dyn(self, X, t, np.concatenate((Fb, Mb)))

    def trim(self, args={}, UengFun=Ueng_of_throttle, UsfcFun=Usfc_of_elevator, debug=False):
        va, gamma, h = args.get('va', self.P.Vref), args.get('gamma', 0), args.get('h', 0)
        if debug:
            print("searching for constant path trajectory with")
            print("  h      {:f} m".format(h))
            print("  va     {:f} m/s".format(va))
            print("  gamma  {:f} deg".format(np.rad2deg(gamma)))
        rho = patm.get_rho(h)
        Pdyn, beta, omega = 0.5*rho*va**2, 0, [0, 0, 0]
        def err_func(args):
            ele, throttle, alpha = args
            Ueng, Usfc = UengFun(throttle, self.P), UsfcFun(ele, self.P)
            eulers = [0, gamma+alpha, 0]
            Fb, f_eng_body = get_forces_body(rho, va, alpha, beta, eulers, omega, Pdyn, Ueng, Usfc, self.P)
            f_weight_body = np.dot(pfw.earth_to_body(eulers), [0., 0., self.P.m*self.P.g])
            Fb = Fb + f_weight_body
            Mb = get_moments_body(alpha, beta, eulers, omega, Pdyn, f_eng_body, Usfc, self.P)
            return [Fb[0], Fb[2], Mb[1]]
        sol = scipy.optimize.root(err_func, [0., 0.5, 0.], method='hybr')
        ele, throttle, alpha = sol.x
        if debug:
            print("""result:
  throttle        : {:f} %
  elevator        : {:f} deg
  angle of attack : {:f} deg""".format(100.*throttle, np.rad2deg(ele), np.rad2deg(alpha)) )

        Ue = np.zeros(self.P.input_nb)
        Ue[:self.P.eng_nb] = UengFun(throttle, self.P)
        Ue[self.P.eng_nb:self.P.input_nb] = UsfcFun(ele, self.P)
        Xe = [0., 0., -h, va, alpha, beta, 0., gamma+alpha, 0., 0., 0., 0.]
        return Xe, Ue


    def plot(self, time, X, U=None, figure=None, window_title="Trajectory"):
        Ueng, Usfc = U[:,:self.P.eng_nb], U[:,self.P.eng_nb:]
        Usolid = np.zeros((len(time), pat_dyn.SolidDM3.iv_size))
        figure = pat_dyn.SolidDM3.plot(self, time, X, None, extra_rows=1)
        ax = plt.subplot(5, 3, 13)
        for i in range(self.P.eng_nb):
            plt.plot(time, 100*Ueng[:,i], label='prop_{}'.format(i))
        ppu.decorate(ax, "Propulsion", xlab="time", ylab="%", legend=True)
        ax = plt.subplot(5, 3, 14)
        for i in range(self.P.sfc_nb):
            plt.plot(time, np.rad2deg(Usfc[:,i]), label='sfc_{}'.format(self.P.sfc_name[i]))
        ppu.decorate(ax, "Surfaces", legend=True, xlab="time", ylab="deg")
        return figure
