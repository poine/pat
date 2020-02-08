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
import pat3.vehicles.fixed_wing.aero as p3_fw_aero


def load(type, variant, param=None):
    mod_str = "pat.vehicles.{:s}.dynamic_model_{:s}".format(type, variant)
    dm = importlib.import_module(mod_str)
    return dm.DynamicModel(param)
#
# TODO
#
# def trim_banked(P, h=0., va=12., throttle=0., phi=0., report=True, debug=False):
#     if report:
#         print("searching for banked constant throttle trajectory with")
#         print("  roll  {:.1f} deg".format(np.rad2deg(phi)))
#         print("  h      {:.1f} m".format(h))
#         print("  va     {:.2f} m/s".format(va))
#         print("  throttle  {:.1f} %".format(100*throttle))
#         R = va**2/np.tan(phi)
#         psi = 0. # does not depend on yaw... or?
#         phid_ref, thetad_ref, psid_ref = 0., 0., va/R
#         #euler_ref, eulerd_ref = [phi, theta, psi], [0, 0, psi_dot]
#         #pref, qref, rref = 0, 0, va/R
#         print(" psid {} deg/s {}s/rev".format(np.rad2deg(psid_ref), 360./np.rad2deg(psid_ref)))

        
#     def err_func(args):
#         gamma, elevator, alpha = args
#         theta = gamma + alpha
#         psi_dot = va/R # not really....
#         euler, eulerd = [phi, theta, psi], [0, 0, psi_dot]
#         #p, q, r = 0, 0, va/R # FIXME...
#         #p, q, r = p3_alg.rvel_of_eulerd(euler, eulerd)
#         p, q, r = 0, np.sin(phi)*np.cos(theta)*va/R, np.cos(phi)*np.cos(theta)*va/R
#         X=[0., 0., -h, va, alpha, 0., phi,  theta, psi, p, q, r]
#         U = np.zeros(P.input_nb)
#         U[0:P.eng_nb] = throttle; U[P.eng_nb+iv_de] = elevator
#         Xdot = dyn(X, 0., U, P)
#         Xdot_ref = [va*math.cos(gamma), 0., -va*math.sin(gamma), 0., 0., 0., phid_ref, thetad_ref, psid_ref, 0, 0, 0]
#         pdb.set_trace()
#         return np.linalg.norm(Xdot - Xdot_ref)

#     p0 = [0., np.deg2rad(-2.), np.deg2rad(1.)]
#     gamma_e, ele_e, alpha_e = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-9)
#     if report:
#         print("""result:
#   gamma           : {:.2f} deg
#   elevator        : {:.2f} deg
#   angle of attack : {:.2f} deg""".format(np.rad2deg(gamma_e), np.rad2deg(ele_e), np.rad2deg(alpha_e)))
#     return gamma_e, ele_e, alpha_e    


import pat3.vehicles.fixed_wing.legacy_parameter
class ParamOld(pat3.vehicles.fixed_wing.legacy_parameter.Param):
    pass
import pat3.vehicles.fixed_wing.simple_6dof_fdm_param
class Param(pat3.vehicles.fixed_wing.simple_6dof_fdm_param.Param):
    pass

#
# Generic
#
class FWDynamicModel:
    def __init__(self, params=None):
        if params is None: params = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
        self.P = Param(params)
        #self.P = ParamOld(params)
        self.X = np.zeros(self.get_sv_size())
        self.X_act = np.zeros(self.input_nb())
        self.T_w2b = np.eye(4)
        self.t, self.dt = 0., 0.01
        self.reset()

    def param(self): return str(self.P)
    def get_sv_size(self): return self.sv_size

    def reset(self, X0=None, t0=0, X_act0=None):
        self.X = np.asarray(X0) if X0 is not None else self.default_state()
        self.X_act = np.asarray(X_act0) if X_act0 is not None else np.zeros(self.input_nb())
        self.t = t0
        self._update_byproducts()
        return self.X

    def run(self, dt, tf, U, atm, act_dyn=False):
        if act_dyn: # actuator dynamics
            act_lambda = -1./np.array([0.1, 0.02, 0.02, 0.02, 0.2])
            self.X_act += dt*((self.X_act-self._clip_input(U))*act_lambda)
        else:
            self.X_act = self._clip_input(U)
       
        foo, self.X = scipy.integrate.odeint(self.dyn, self.X, [self.t, self.t+dt], args=(self.X_act, atm))#, hmax=0.001)
        self.t += dt
        self._update_byproducts()
        return self.X
    

    #
    # Actuator allocation, still needs some love
    #
    #_iv_th   = 0    # throttle
    _iv_da   = 0    # aileron
    _iv_de   = 1    # elevator
    _iv_dr   = 2    # rudder
    _iv_df   = 3    # flap

    def input_nb(self): return self.P.input_nb

    def _clip_input(self, U):
        Uclip = np.copy(U)
        # clip input
        max_ail, max_ele, max_rud = np.deg2rad(30), np.deg2rad(25), np.deg2rad(30)
        Uclip[self.iv_dth()] = np.clip(U[self.iv_dth()], 0., 1.)
        Uclip[self.iv_da()]  = np.clip(U[self.iv_da()], -max_ail, max_ail)
        Uclip[self.iv_de()]  = np.clip(U[self.iv_de()], -max_ele, max_ele)
        Uclip[self.iv_dr()]  = np.clip(U[self.iv_dr()], -max_rud, max_rud)
        return Uclip

    def iv_dth(self):
        if self.P.eng_nb>1: return range(0,self.P.eng_nb)
        else: return 0
    def iv_da(self): return self.P.eng_nb + self._iv_da
    def iv_de(self): return self.P.eng_nb + self._iv_de
    def iv_dr(self): return self.P.eng_nb + self._iv_dr
    def iv_df(self): return self.P.eng_nb + self._iv_df


    #
    # Trims (using forces and moments rather than full dynamic equation)
    #
    def trim(self, args, report=True, debug=False, aero=False, atm=None):
        '''
        Find a set of state and input variables corresponding
        to another set of state and input variables specified as args
        (yes, really...)
        '''
        if report: print('### FWDynamicModel')
        flaps = args.get('flaps', 0.)
        phi = 0.
        ail_e, rud_e = 0., 0.
        if 'va' in args and 'gamma' in args :
            h, va, gamma_e = args['h'], args['va'], args['gamma']
            if aero:
                thr_e, ele_e, alpha_e = self.trim_aero_cst_path(h, va, gamma_e, flaps, report, debug)
            else:
                thr_e, ele_e, alpha_e = self.trim_cst_path(h, va, gamma_e, flaps, report, debug)
        elif 'va' in args and 'throttle' in args and 'phi' in args:
            h, va, thr_e, phi = args['h'], args['va'], args['throttle'], args['phi']
            gamma_e, ele_e, ail_e, rud_e, alpha_e = self.trim_aero_banked_cst_throttle(h, va, thr_e, flaps, phi, report, debug, atm)
        elif 'va' in args and 'throttle' in args:
            h, va, thr_e = args['h'], args['va'], args['throttle']
            if aero:
                gamma_e, ele_e, alpha_e = self.trim_aero_cst_throttle(h, va, thr_e, flaps, report, debug)
            else:
                gamma_e, ele_e, alpha_e = self.trim_cst_throttle(h, va, thr_e, flaps, report, debug)
        else:
            print('trim not possible. Specify h, and either va and gamma, or va and throttle')
            
        Ue = np.zeros(self.input_nb())
        Ue[self.iv_dth()] = thr_e; Ue[self.iv_de()] = ele_e; Ue[self.iv_df()] = flaps
        Ue[self.iv_da()] = ail_e; Ue[self.iv_dr()] = rud_e; 
        theta_e = gamma_e+alpha_e
        psid = self.P.g*np.tan(phi)/va
        p_e, q_e, r_e = 0., np.sin(phi)*np.cos(theta_e)*psid, np.cos(phi)*np.cos(theta_e)*psid
        Xe = [0., 0., -h, va, alpha_e, 0., phi, theta_e, 0., p_e, q_e, r_e]
        _Xe = self.from_six_dof_aero_euler(Xe)
        return _Xe, Ue

    def trim_aero_cst_path(self, h=0., va=12., gamma=0., flaps=0., report=True, debug=False):
        if report:
            print("searching (aero) for constant path trajectory with")
            print("  h      {:.1f} m".format(h))
            print("  va     {:.2f} m/s".format(va))
            print("  gamma  {:.2f} deg".format(np.rad2deg(gamma)))
            print("  flaps  {:.2f} deg".format(np.rad2deg(flaps)))
        beta = 0.
        rho = p3_atm.get_rho(h)
        Pdyn = 0.5*rho*va**2
        X_rvel_body = [0, 0, 0]
        def err_func(args):
            throttle, elevator, alpha = args
            Ueng, Usfc = [throttle], [0, elevator, 0, flaps]
            f_aero_body = p3_fw_aero.get_f_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
            f_eng_body  =  p3_fw_aero.get_f_eng_body(h, va, Ueng, self.P)

            eul = phi, theta, psi = [0, gamma+alpha, 0]
            earth_to_body_R = p3_alg.rmat_of_euler(eul)
            f_weight_body = np.dot(earth_to_body_R, [0., 0., self.P.m*self.P.g])
            forces_body = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body
            m_aero_body = p3_fw_aero.get_m_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
            m_eng_body = p3_fw_aero.get_m_eng_body(f_eng_body, self.P)
            m_body = m_aero_body + m_eng_body
            return np.linalg.norm(np.concatenate((forces_body, m_body)))
        
        p0 = [0.2, np.deg2rad(2.), np.deg2rad(0.)]
        thr_e, ele_e, alpha_e = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-9)
        if report:
            print("""result:
  throttle        : {:.2f} %
  elevator        : {:.2f} deg
  angle of attack : {:.2f} deg""".format(100.*thr_e, np.rad2deg(ele_e), np.rad2deg(alpha_e)))
        return thr_e, ele_e, alpha_e    

    def trim_aero_cst_throttle(self, h=0., va=12., throttle=0., flaps=0., report=True, debug=False, atm=None):
        if report:
            print("searching (aero) for constant throttle trajectory with")
            print("  h      {:.1f} m".format(h))
            print("  va     {:.2f} m/s".format(va))
            print("  throttle  {:.1f} %".format(100*throttle))
            print("  flap  {:.1f} deg".format(np.rad2deg(flaps)))
        beta = 0.; rho = p3_atm.get_rho(h); Pdyn = 0.5*rho*va**2
        X_rvel_body = [0, 0, 0]
        def err_func(args):
            gamma, elevator, alpha = args
            Ueng, Usfc = [throttle], [0, elevator, 0, flaps]
            f_aero_body = p3_fw_aero.get_f_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
            f_eng_body  =  p3_fw_aero.get_f_eng_body(h, va, Ueng, self.P)
            eul = phi, theta, psi = [0, gamma+alpha, 0]
            earth_to_body_R = p3_alg.rmat_of_euler(eul)
            f_weight_body = np.dot(earth_to_body_R, [0., 0., self.P.m*self.P.g])
            forces_body = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body
            m_aero_body = p3_fw_aero.get_m_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
            m_eng_body = p3_fw_aero.get_m_eng_body(f_eng_body, self.P)
            m_body = m_aero_body + m_eng_body
            return np.linalg.norm(np.concatenate((forces_body, m_body)))
        
        p0 = [0., np.deg2rad(-2.), np.deg2rad(1.)]
        gamma_e, ele_e, alpha_e = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-9)
        if report:
            print("""result:
  gamma           : {:.2f} deg
  elevator        : {:.2f} deg
  angle of attack : {:.2f} deg""".format(np.rad2deg(gamma_e), np.rad2deg(ele_e), np.rad2deg(alpha_e)))
        return gamma_e, ele_e, alpha_e    


    # http://zeta.math.utsa.edu/~gokhman/ftp/courses/notes/frenet.pdf
    def trim_aero_banked_cst_throttle(self, h=0., va=12., throttle=0., flaps=0., phi=0., report=True, debug=False, atm=None):
        beta = 0.; rho = p3_atm.get_rho(h); Pdyn = 0.5*rho*va**2
        psi, psid = 0., self.P.g*np.tan(phi)/va
        def err_func(args):
            theta, elevator, aileron, rudder, alpha = args
            Ueng, Usfc = [throttle], [aileron, elevator, rudder, flaps]
            X_rvel_body = p, q, r = 0., np.sin(phi)*np.cos(theta)*psid, np.cos(phi)*np.cos(theta)*psid
            #X_rvel_body = p, q, r = 0, np.sin(phi)*np.cos(theta)*va/R, np.cos(phi)*np.cos(theta)*va/R
            wind_ned = [0, 0, 0]
            earth_to_body_R = p3_alg.rmat_of_euler([phi, theta, psi])
            ivel_b = p3_fr.vel_aero_to_body_R([va, alpha, beta], earth_to_body_R, wind_ned)
            f_aero_body = p3_fw_aero.get_f_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
            f_eng_body  =  p3_fw_aero.get_f_eng_body(h, va, Ueng, self.P)
            #f_body  = f_aero_body + np.sum(f_eng_body, axis=0)
            #f_ned = np.dot(earth_to_body_R.T, f_body).squeeze()
            #f_ned[2] += self.P.m*self.P.g
            #f_ned[1] -= self.P.m*va*psid 
            f_weight_body = np.dot(earth_to_body_R, [0., 0., self.P.m*self.P.g])
            f_body = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body
            uvwd = 1/self.P.m*f_body - np.cross(X_rvel_body, ivel_b)

            m_aero_body = p3_fw_aero.get_m_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
            m_eng_body = p3_fw_aero.get_m_eng_body(f_eng_body, self.P)
            m_body = m_aero_body + m_eng_body
            #pdb.set_trace()
            #pqrd = np.dot(self.P.invI, m_body - np.cross(X_rvel_body, np.dot(self.P.I, X_rvel_body)))
            pqrd = m_body - np.cross(X_rvel_body, np.dot(self.P.I, X_rvel_body))
            foo = np.concatenate((uvwd, pqrd))
            #return np.linalg.norm(np.concatenate((f_ned, m_body)))
            #print foo
            return np.linalg.norm(foo)#+100.*abs(aileron)
            
        p0 = [0., np.deg2rad(-2.), np.deg2rad(-1.), np.deg2rad(2.), 0.]
        theta, elevator, aileron, rudder, alpha = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-12)
        gamma = theta-alpha
        if report:
            print("""result:
  gamma           : {:.2f} deg
  theta           : {:.2f} deg
  elevator        : {:.2f} deg
  aileron         : {:.2f} deg
  rudder          : {:.2f} deg
  angle of attack : {:.2f} deg""".format(np.rad2deg(gamma), np.rad2deg(theta), np.rad2deg(elevator), np.rad2deg(aileron), np.rad2deg(rudder), np.rad2deg(alpha))) 
        return gamma, elevator, aileron, rudder, alpha 

    def foo_trim_aero_banked_cst_throttle(self, h=0., va=12., throttle=0., flaps=0., phi=0., report=True, debug=False, atm=None):
        if report:
            print("searching (aero) for banked/constant throttle trajectory with")
            print("  h      {:.1f} m".format(h))
            print("  va     {:.2f} m/s".format(va))
            print("  phi    {:.1f} deg".format(np.rad2deg(phi)))
            print("  throttle  {:.1f} %".format(100*throttle))
            print("  flap  {:.1f} deg".format(np.rad2deg(flaps)))
        gamma, elevator, aileron, rudder, alpha = self.trim_aero_banked_cst_throttle(h, va, throttle, flaps, phi, report, debug, atm)
        Ue = np.zeros(self.input_nb())
        Ue[self.iv_dth()] = throttle; Ue[self.iv_de()] = elevator; Ue[self.iv_df()] = flaps
        Ue[self.iv_dr()] = rudder; Ue[self.iv_da()] = aileron
        psid = self.P.g*np.tan(phi)/va
        theta = gamma+alpha
        p, q, r = 0., np.sin(phi)*np.cos(theta)*psid, np.cos(phi)*np.cos(theta)*psid
        Xe = [0., 0., -h, va, alpha, 0., phi, theta, 0.,  p, q, r]
        _Xe = self.from_six_dof_aero_euler(Xe)

        #_Xe_be = p3_fr.SixDOFBodyEuler.state_from_ee(_Xe)
        #print _Xe_be
        #print p3_fr.SixDOFBodyEuler.
        return _Xe, Ue

#
#
# Aero/Euler
#
#    
class DynamicModel(p3_fr.SixDOFAeroEuler, FWDynamicModel):

     
    def __init__(self, params=None):
        print("Info: Dynamic fixed wing legacy")
        FWDynamicModel.__init__(self, params)

    def dyn(self, X, t, U, atm=None):
        Ueng = U[self.P.u_slice_eng()]                    # engines part of input vector
        Usfc = U[self.P.u_slice_sfc()]                    # control surfaces part of input vector
        X_pos = X[self.sv_slice_pos]                      # ned pos
        X_avel = va, alpha, beta = X[self.sv_slice_vaero] # airvel, alpha, beta
        X_euler = X[self.sv_slice_eul]                    # euler angles                          
        X_rvel_body = X[self.sv_slice_rvel]               # body rotational velocities
        earth_to_body_R = p3_alg.rmat_of_euler(X_euler)
        if 0:
            R_aero_to_body = p3_fr.R_aero_to_body(alpha, beta)
            wind_ned = atm.get_wind(X_pos, t) if atm is not None else [0, 0, 0]
            waccel_body = [0, 0, 0]  # np.dot(earth_to_body, waccel_earth)
            avel_aero = [X[sv_v], 0., 0.]
            avel_body = np.dot(R_aero_to_body, avel_aero)
            ivel_world = p3_fr.vel_aero_to_world_euler(X_avel, X_euler, wind_ned)
            ivel_body = np.dot(earth_to_body_R, ivel_world)
        
        rho = p3_atm.get_rho(-X[self.sv_z])
        Pdyn = 0.5*rho*va**2

        # Forces
        f_aero_body = p3_fw_aero.get_f_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
        f_eng_body  =  p3_fw_aero.get_f_eng_body(-X[self.sv_z], va, Ueng, self.P)
        f_weight_body = np.dot(earth_to_body_R, [0., 0., self.P.m*self.P.g])
        forces_body = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body
        # Moments
        m_aero_body = p3_fw_aero.get_m_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
        m_eng_body = p3_fw_aero.get_m_eng_body(f_eng_body, self.P)
        m_body = m_aero_body + m_eng_body
        
        if 1:
            return p3_fr.SixDOFAeroEuler.cont_dyn(X, t, np.concatenate((forces_body, m_body)), self.P, atm)
        else:
            Xdot = np.zeros(p3_fr.SixDOFAeroEuler.sv_size)
            # Translational kinematics
            Xdot[sv_x:sv_z+1] = ivel_world
            # Translational dynamics
            iaccel_body = 1./self.P.m*forces_body - np.cross(X_rvel_body, ivel_body)
            aaccel_body = iaccel_body - waccel_body
            Xdot[sv_v] = np.inner(avel_body, aaccel_body)/X[self.sv_va]
            (avx, avy, avz), (aax, aay, aaz) = avel_body, aaccel_body
            Xdot[sv_alpha] = (avx*aaz - avz*aax)/(avx**2+avz**2)
            Xdot[sv_beta] = (X[self.sv_va]*aay - avy*Xdot[self.sv_va]) / X[self.sv_va] / math.sqrt(avx**2+aaz**2)
            # Rotational kinematics
            Xdot[self.sv_slice_eul] = p3_alg.euler_derivatives(X_euler, X_rvel_body)
            # Rotational dynamics
            raccel_body = np.dot(self.P.invI, m_body - np.cross(X_rvel_body, np.dot(self.P.I, X_rvel_body)))
            Xdot[self.sv_slice_rvel] = raccel_body
            return Xdot
            

    def _update_byproducts(self):
        self.T_w2b[:3,:3] = p3_alg.rmat_of_euler(self.X[self.sv_slice_eul]).T
        self.T_w2b[:3,3] = self.X[self.sv_slice_pos]

    #def get_sv_size(self): return p3_fr.SixDOFAeroEuler.sv_size
 
    def name(self): return "Fixed Wing, Aero/Euler ({:s})".format(self.P.name)

    def default_state(self): return np.array([0, 0, 0,   10, 0, 0,   0, 0, 0,   0, 0, 0])

    # def thrust_of_throttle(self, throttle, i, h, v):
    #     rho = p3_atm.get_rho(h)
    #     return thrust_of_throttle(throttle, self.P.fmaxs[i], 
    #                               rho, self.P.rhois[i], self.P.nrhos[i], 
    #                               v, self.P.Vis[i], self.P.nVs[i])

    def state_as_six_dof_euclidian_euler(self, X, atm=None, t=0.):
        return p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(X, atm, t)
    
    def get_jacobian(self, Xe, Ue):
        A,B = p3_u.num_jacobian(Xe, Ue, None, self.dyn)
        return A, B 

    def trim_cst_path(self, h=0., va=12., gamma=0., flaps=0., report=True, debug=False, atm=None):
        if report:
            print("searching (dyn) for constant path trajectory with")
            print("  h      {:.1f} m".format(h))
            print("  va     {:.2f} m/s".format(va))
            print("  gamma  {:.2f} deg".format(np.rad2deg(gamma)))
            print("  flaps  {:.2f} deg".format(np.rad2deg(flaps)))

        def err_func(args):
            throttle, elevator, alpha = args
            X=[0., 0., -h, va, alpha, 0., 0.,  gamma+alpha, 0., 0., 0., 0.]
            U = np.zeros(self.input_nb())
            U[self.iv_dth()] = throttle; U[self.iv_de()] = elevator; U[self.iv_df()] = flaps
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

    def trim_cst_throttle(self, h=0., va=12., throttle=0., flaps=0., report=True, debug=False, atm=None):
        if report:
            print("searching (dyn) for constant throttle trajectory with")
            print("  h      {:.1f} m".format(h))
            print("  va     {:.2f} m/s".format(va))
            print("  throttle  {:.1f} %".format(100*throttle))
            print("  flap  {:.1f} deg".format(np.rad2deg(flaps)))
        def err_func(args):
            gamma, elevator, alpha = args
            X=[0., 0., -h, va, alpha, 0., 0.,  gamma+alpha, 0., 0., 0., 0.]
            U = np.zeros(self.input_nb())
            U[self.iv_dth()] = throttle; U[self.iv_de()] = elevator; U[self.iv_df()] = flaps
            Xdot = self.dyn(X, 0., U, atm)
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
        


    def state_str(self, X=None):
        return p3_fr.SixDOFAeroEuler.state_str(X if X is not None else self.X)

    def plot_trajectory(self, time, Xae, U=None, figure=None, window_title="Trajectory", legend=None, label='', filename=None, atm=None): 
        plot_trajectory_ae(time, Xae, U, figure, window_title, legend, label, filename, atm)

    def plot_trajectory_as_ae(self, time, X, U=None, figure=None, window_title="Trajectory", legend=None, label='', filename=None, atm=None): 
        plot_trajectory_ae(time, X, U, figure, window_title, legend, label, filename, atm)

    def plot_trajectory_as_ee(self, time, Xae, U=None, figure=None, window_title="Trajectory", legend=None, label='', filename=None, atm=None): 
        Xee = np.array([p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(_X, atm, _t) for _X, _t in zip(Xae, time)])
        plot_trajectory_ee(time, Xee, U, figure, window_title, legend, label, filename, atm)


        
#
#
# Euclidian/Euler
#
#    
class DynamicModel_ee(p3_fr.SixDOFEuclidianEuler, FWDynamicModel):

    def __init__(self, params_filesname=None):
        print("Info: Dynamic fixed wing ee")
        FWDynamicModel.__init__(self, params_filesname)
      

    def dyn(self, X, t, U, atm):

        Ueng = U[self.P.u_slice_eng()]                    # engines part of input vector
        Usfc = U[self.P.u_slice_sfc()]                    # control surfaces part of input vector
        X_pos = X[self.sv_slice_pos]                      #
        X_vel = X[self.sv_slice_vel]                      #
        X_euler = X[self.sv_slice_eul]                    # euler angles                          
        X_rvel_body = X[self.sv_slice_rvel]               # body rotational velocities

        earth_to_body_R = p3_alg.rmat_of_euler(X_euler)
        wind_ned = atm.get_wind(X_pos, t) if atm is not None else [0, 0, 0]
        va, alpha, beta = p3_fr.vel_world_to_aero_eul(X_vel, X_euler, wind_ned)
        h = -X[self.sv_z]
        rho = p3_atm.get_rho(h)
        Pdyn = 0.5*rho*va**2

        f_aero_body = p3_fw_aero.get_f_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
        f_eng_body  =  p3_fw_aero.get_f_eng_body(h, va, Ueng, self.P)
        f_weight_body = np.dot(earth_to_body_R, [0., 0., self.P.m*self.P.g])
        forces_body = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body

        m_aero_body = p3_fw_aero.get_m_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
        m_eng_body = p3_fw_aero.get_m_eng_body(f_eng_body, self.P)
        moments_body = m_aero_body + m_eng_body
        
        if 1:
            _U = np.concatenate((forces_body, moments_body))
            return p3_fr.SixDOFEuclidianEuler.cont_dyn(X, t, _U, self.P)
        else:
            Xdot = np.zeros(self.sv_size)
            # Position kinematics
            Xdot[self.sv_slice_pos] = X[self.sv_slice_vel]
            # Newton for forces in world frame
            forces_world = np.dot(earth_to_body_R.T, forces_body)
            Xdot[self.sv_slice_vel] = 1./self.P.m * forces_world
            # Orientation kinematics
            Xdot[self.sv_slice_eul] = p3_alg.euler_derivatives(X_euler, X_rvel_body)
            # Newton for moments in body frame
            raccel_body = np.dot(self.P.invI, moments_body - np.cross(X_rvel_body, np.dot(self.P.I, X_rvel_body)))
            Xdot[self.sv_slice_rvel] = raccel_body
            return Xdot

    def _update_byproducts(self):
        self.T_w2b[:3,:3] = p3_alg.rmat_of_euler(self.X[self.sv_slice_eul]).T
        self.T_w2b[:3,3] = self.X[self.sv_slice_pos]
        
    def default_state(self): return np.array([0, 0, 0,   10, 0, 0,   0, 0, 0,   0, 0, 0])
    
    def trim_cst_path(self, h, va, gamma, flaps, report, debug, atm=None):
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

    def trim_cst_throttle(self, h, va, throttle, flaps, report, debug, atm=None):
        if report:
            print("searching for constant throttle trajectory with")
            print("  h      {:.1f} m".format(h))
            print("  va     {:.2f} m/s".format(va))
            print("  throttle  {:.2f} %".format(100.*throttle))
        wvel_world = (atm.get_wind([0, 0, -h], t=0) if atm is not None else [0, 0, 0])
        def err_func(args):
            gamma, elevator, alpha = args
            xd, yd, zd = va*math.cos(gamma), 0., -va*math.sin(gamma)
            phi, theta, psi = 0.,  gamma+alpha, 0.
            X=[0., 0., -h, xd, yd, zd, phi, theta, psi, 0., 0., 0.]
            U = np.zeros(self.P.input_nb)
            U[self.P.u_slice_eng()], U[self.P.u_elevator()] = throttle, elevator
            Xdot = self.dyn(X, 0., U, atm)
            Xdot_ref = [va*math.cos(gamma), 0., -va*math.sin(gamma), 0., 0., 0., 0., 0., 0., 0., 0., 0.]
            return np.linalg.norm(Xdot - Xdot_ref)

        p0 = [np.deg2rad(0.), np.deg2rad(2.), np.deg2rad(0.)]
        gamma_e, ele_e, alpha_e = scipy.optimize.fmin_powell(err_func, p0, disp=debug, ftol=1e-9)

        if report:
            print("""result:
    gamma           : {:.2f} deg
    elevator        : {:.2f} deg
    angle of attack : {:.2f} deg""".format(np.rad2deg(gamma_e), np.rad2deg(ele_e), np.rad2deg(alpha_e)))

        return gamma_e, ele_e, alpha_e

    #def get_sv_size(self): return self.sv_size
 
    def plot_trajectory(self, time, Xee, U=None, figure=None, window_title="Trajectory", legend=None, label='', filename=None, atm=None): 
        return plot_trajectory_ee(time, Xee, U, figure, window_title, legend, label, filename, atm)

    def plot_trajectory_as_ee(self, time, Xee, U=None, figure=None, window_title="Trajectory", legend=None, label='', filename=None, atm=None): 
        return plot_trajectory_ee(time, Xee, U, figure, window_title, legend, label, filename, atm)

    def plot_trajectory_as_ae(self, time, Xae, U=None, figure=None, window_title="Trajectory", legend=None, label='', filename=None, atm=None):
        Xae = np.array([p3_fr.SixDOFEuclidianEuler.to_six_dof_aero_euler(_X, atm, _t) for _X, _t in zip(Xae, time)])
        return plot_trajectory_ae(time, Xae, U, figure, window_title, legend, label, filename, atm)

    def state_as_six_dof_euclidian_euler(self, X, atm=None, t=0.):
        return X
    def state_as_six_dof_aero_euler(self, X, atm=None, t=0.):
        return p3_fr.SixDOFEuclidianEuler.to_six_dof_aero_euler(X, atm, t)
    # FIXME
    def get_jacobian(self, Xe, Ue):
        A,B = p3_u.num_jacobian(Xe, Ue, None, self.dyn)
        return A, B 


#
#
# Euclidian/Quaternion
#
#    
class DynamicModel_eq(p3_fr.SixDOFEuclidianQuat, FWDynamicModel):

    def __init__(self, params_filesname=None):
        print("Info: Dynamic fixed wing eq")
        FWDynamicModel.__init__(self, params_filesname)
        
    def default_state(self): return np.array([0, 0, 0,   10, 0, 0,   1, 0, 0, 0,   0, 0, 0])
    #def get_sv_size(self): return self.sv_size
    def _update_byproducts(self):
        self.T_w2b[:3,:3] = p3_alg.rmat_of_quat(self.X[self.sv_slice_quat]).T
        self.T_w2b[:3,3] = self.X[self.sv_slice_pos]
    def trim_cst_path(self, h, va, gamma, flaps, report, debug, atm=None):
        thr_e, ele_e, alpha_e = 0.17, -0.05, 0.01
        return thr_e, ele_e, alpha_e
    def state_as_six_dof_euclidian_euler(self, X, atm=None):
        return p3_fr.SixDOFEuclidianQuat.to_six_dof_euclidian_euler(X)
    def state_as_six_dof_aero_euler(self, X, atm=None):
        return p3_fr.SixDOFEuclidianQuat.to_six_dof_aero_euler(X)

    def dyn(self, X, t, U, atm):
        Ueng = U[self.P.u_slice_eng()]                    # engines part of input vector
        Usfc = U[self.P.u_slice_sfc()]                    # control surfaces part of input vector
        X_pos = X[self.sv_slice_pos]                      #
        X_vel = X[self.sv_slice_vel]                      #
        X_quat = X[self.sv_slice_quat]                    # quaternion                          
        X_rvel_body = X[self.sv_slice_rvel]               # body rotational velocities

        earth_to_body_R = p3_alg.rmat_of_quat(X_quat)
        wind_ned = atm.get_wind(X_pos, t) if atm is not None else [0, 0, 0]
        va, alpha, beta = p3_fr.vel_world_to_aero_quat(X_vel, X_quat, wind_ned)
        h = -X[self.sv_z]; rho = p3_atm.get_rho(h); Pdyn = 0.5*rho*va**2

        f_aero_body = p3_fw_aero.get_f_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
        f_eng_body  =  p3_fw_aero.get_f_eng_body(h, va, Ueng, self.P)
        f_weight_body = np.dot(earth_to_body_R, [0., 0., self.P.m*self.P.g])
        forces_body = f_aero_body + np.sum(f_eng_body, axis=0) + f_weight_body

        m_aero_body = p3_fw_aero.get_m_aero_body(va, alpha, beta, X_rvel_body, Usfc, self.P, Pdyn)
        m_eng_body = p3_fw_aero.get_m_eng_body(f_eng_body, self.P)
        moments_body = m_aero_body + m_eng_body
        
        _U = np.concatenate((forces_body, moments_body))
        return p3_fr.SixDOFEuclidianQuat.cont_dyn(X, t, _U, self.P)

    def plot_trajectory(self, time, Xeq, U=None, figure=None, window_title="Trajectory", legend=None, label='', filename=None, atm=None): 
        return self.plot_trajectory_as_ee(time, Xeq, U, figure, window_title, legend, label, filename, atm)

    def plot_trajectory_as_ee(self, time, Xeq, U=None, figure=None, window_title="Trajectory", legend=None, label='', filename=None, atm=None):
        Xee = np.array([p3_fr.SixDOFEuclidianQuat.to_six_dof_euclidian_euler(_X, atm, _t) for _X, _t in zip(Xeq, time)])
        plot_trajectory_ee(time, Xee, U, figure, window_title, legend, label, filename, atm)
        
#
# Some plotting functions
#
def plot_trajectory_ae(time, X, U=None, figure=None, window_title="Trajectory",
                       legend=None, label='', filename=None, atm=None):
    s = p3_fr.SixDOFAeroEuler
    margins=(0.04, 0.05, 0.98, 0.96, 0.20, 0.34)
    figure = p3_pu.prepare_fig(figure, window_title, figsize=(20.48, 10.24), margins=margins)

    plots = [("x", "m", X[:,s.sv_x]), ("y", "m", X[:,s.sv_y]), ("z", "m", X[:,s.sv_z]),
             ("v",         "m/s", X[:,s.sv_va]),
             ("$\\alpha$", "deg", np.rad2deg(X[:,s.sv_alpha])),
             ("$\\beta$",  "deg", np.rad2deg(X[:,s.sv_beta])),
             ("$\phi$",    "deg", np.rad2deg(X[:,s.sv_phi])),
             ("$\\theta$", "deg", np.rad2deg(X[:,s.sv_theta])),
             ("$\\psi$",   "deg", np.rad2deg(X[:,s.sv_psi])),
             ("$p$",     "deg/s", np.rad2deg(X[:,s.sv_p])),
             ("$q$",     "deg/s", np.rad2deg(X[:,s.sv_q])),
             ("$r$",     "deg/s", np.rad2deg(X[:,s.sv_r]))]
    
    nrow = 5 if U is not None else 4
    for i, (title, ylab, data) in enumerate(plots):
        ax = plt.subplot(nrow, 3, i+1)
        plt.plot(time, data, label=label)
        p3_pu.decorate(ax, title=title, ylab=ylab)
            
    if legend!=None:
        plt.legend(legend, loc='best')

    for i, min_yspan in enumerate([1., 1., 1.,  1., 1., 1.,  1., 1., 1.,  1., 1., 1.]):
        p3_pu.ensure_yspan(plt.subplot(nrow,3,i+1), min_yspan)

    iv_da, iv_de, iv_dr, iv_df = 1, 2, 3, 4 # FIXME... needs P?
    if U is not None:
        ax = figure.add_subplot(5, 3, 13)
        ax.plot(time, 100*U[:, 0])
        p3_pu.decorate(ax, title="$d_{th}$", ylab="%", min_yspan=1.)
        ax = figure.add_subplot(5, 3, 14)
        ax.plot(time, np.rad2deg(U[:, iv_da]), label="aileron")
        ax.plot(time, np.rad2deg(U[:, iv_dr]), label="rudder")
        p3_pu.decorate(ax, title="$d_a/d_r$", ylab="deg", min_yspan=1., legend=True)
        ax = figure.add_subplot(5, 3, 15)
        ax.plot(time, np.rad2deg(U[:, iv_de]), label="elevator")
        ax.plot(time, np.rad2deg(U[:, iv_df]), label="flap")
        p3_pu.decorate(ax, title="$d_e/d_f$", ylab="deg", min_yspan=1., legend=True)
        
    return figure
        
# FIXME, move that elsewhere
def plot_trajectory_ee(time, Xee, U=None, figure=None, window_title="Trajectory",
                       legend=None, label='', filename=None, atm=None):
    s = p3_fr.SixDOFEuclidianEuler
    margins=(0.04, 0.05, 0.98, 0.96, 0.20, 0.34)
    figure = p3_pu.prepare_fig(figure, window_title, figsize=(20.48, 10.24), margins=margins)
    plots = [("x", "m", Xee[:,s.sv_x]), ("y", "m", Xee[:,s.sv_y]), ("z", "m", Xee[:,s.sv_z]),
             ("$\dot{x}$",  "m/s", Xee[:, s.sv_xd]),
             ("$\dot{y}$",  "m/s", Xee[:, s.sv_yd]),
             ("$\dot{z}$",  "m/s", Xee[:, s.sv_zd]),
             ("$\phi$",     "deg",   np.rad2deg(Xee[:,s.sv_phi])),
             ("$\\theta$",  "deg",   np.rad2deg(Xee[:,s.sv_theta])),
             ("$\\psi$",    "deg",   np.rad2deg(Xee[:,s.sv_psi])),
             ("$p$",        "deg/s", np.rad2deg(Xee[:,s.sv_p])),
             ("$q$",        "deg/s", np.rad2deg(Xee[:,s.sv_q])),
             ("$r$",        "deg/s", np.rad2deg(Xee[:,s.sv_r]))]

    nrow = 5 if U is not None else 4
    for i, (title, ylab, data) in enumerate(plots):
        ax = plt.subplot(nrow, 3, i+1)
        plt.plot(time, data, label=label)
        p3_pu.decorate(ax, title=title, ylab=ylab, legend=True)
            
    #if legend!=None:
    #    plt.legend(legend, loc='best')

    for i, min_yspan in enumerate([1., 1., 1.,  1., 1., 1.,  1., 1., 1.,  1., 1., 1.]):
        p3_pu.ensure_yspan(plt.subplot(nrow,3,i+1), min_yspan)

    iv_da, iv_de, iv_dr, iv_df = 1, 2, 3, 4 # FIXME... needs P?
    if U is not None:
        ax = figure.add_subplot(5, 3, 13)
        ax.plot(time, 100*U[:, 0])
        p3_pu.decorate(ax, title="$d_{th}$", ylab="%", min_yspan=0.01)
        ax = figure.add_subplot(5, 3, 14)
        ax.plot(time, np.rad2deg(U[:, iv_da]), label="aileron")
        ax.plot(time, np.rad2deg(U[:, iv_dr]), label="rudder")
        p3_pu.decorate(ax, title="$d_a/d_r$", ylab="deg", min_yspan=1., legend=True)
        ax = figure.add_subplot(5, 3, 15)
        ax.plot(time, np.rad2deg(U[:, iv_de]), label="elevator")
        ax.plot(time, np.rad2deg(U[:, iv_df]), label="flap")
        p3_pu.decorate(ax, title="$d_e/d_f$", ylab="deg", min_yspan=1., legend=True)
    
    return figure
