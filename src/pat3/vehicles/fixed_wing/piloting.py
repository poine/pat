import os, sys, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import control.matlab

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.frames as p3_fr
import pat3.atmosphere as p3_atm
import pat3.algebra as p3_alg

#
# Composite roll/pitch reference model
#
class EulerRef:
    GRAVITY_MSS = 9.80665
    
    def __init__(self):
        self.phi_ref =   p3_u.SecOrdLinRef(omega=9., xi=0.9, sats=[np.deg2rad(100.), np.deg2rad(400.)]) # rvel, raccel
        self.theta_ref = p3_u.SecOrdLinRef(omega=5., xi=0.9, sats=[np.deg2rad(45.), np.deg2rad(100.)])  # rvel, raccel
        self.p, self.q, self.r = 0, 0, 0
        
    def reset(self, phi, theta, phi_d=0., theta_d=0., phi_dd=0, theta_dd=0, va=12.):
        psid = self.GRAVITY_MSS*np.tan(phi)/va
        psidd = self.GRAVITY_MSS*phi_d/np.cos(phi)**2/va
        #r_coord = np.cos(phi)*np.cos(theta)*psid
        #phi_dot, theta_dot, psi_dot = p3_alg.euler_derivatives([phi, theta, 0.], [p, q, r_coord])
        self.phi_ref.reset(np.array([phi, phi_d, phi_dd]))
        self.theta_ref.reset(np.array([theta, theta_d, theta_dd]))
        self.euler = np.array([phi, theta, 0])
        
        self.euler_d = np.array([phi_d, theta_d, psid])
        self.euler_dd = np.array([phi_dd, theta_dd, psidd])

        self.rvel_body = p3_alg.rvel_of_eulerd(self.euler, self.euler_d)
        self.raccel_body = p3_alg.raccel_of_eulerdd(self.euler, self.euler_d, self.euler_dd)
        #self.rvel_body = np.array([0, 0, 0]) # FIXME
        #self.raccel_body = np.array([0, 0, 0]) # FIXME
        
    def run(self, dt, phi_sp, theta_sp, va):
        phi_ref, phid_ref, phidd_ref = self.phi_ref.run(dt, phi_sp)
        theta_ref, thetad_ref, thetadd_ref = self.theta_ref.run(dt, theta_sp)
        psid = self.GRAVITY_MSS*np.tan(phi_ref)/va
        psidd = self.GRAVITY_MSS*phid_ref/np.cos(phi_ref)**2/va
        #R = va**2/(g*np.tan(phi))
        #p, q, r = 0, np.sin(phi)*np.cos(theta)*va/R, np.cos(phi)*np.cos(theta)*va/R
        self.euler = np.array([phi_ref, theta_ref, 0])
        self.euler_d = np.array([phid_ref, thetad_ref, psid])
        self.euler_dd = np.array([phidd_ref, thetadd_ref, psidd])
        self.phi = phi_ref
        self.theta = theta_ref
        if 1:
            self.p, self.q, self.r = self.rvel_body = p3_alg.rvel_of_eulerd(self.euler, self.euler_d)
            self.pd, self.qd, self.rd = self.raccel_body = p3_alg.raccel_of_eulerdd(self.euler, self.euler_d, self.euler_dd)
        else:
            self.p = phid_ref 
            #self.q = np.sin(phi_ref)**2/np.cos(phi_ref)*np.cos(theta_ref)*self.GRAVITY_MSS/va
            self.q = np.cos(phi_ref)*thetad_ref + np.sin(phi_ref)*np.cos(theta_ref)*psid
            #self.r = np.sin(phi_ref)*np.cos(theta_ref)*self.GRAVITY_MSS/va
            self.r = -np.sin(phi_ref)*thetad_ref + np.cos(phi_ref)*np.cos(theta_ref)*psid

class EulerRefLogger():
    def __init__(self):
        self.eulers = []
        self.euler_ds = []
        self.euler_dds = []
        self.rvel_bodys = []
        self.raccel_bodys = []
        
    def log(self, _ref):
        self.eulers.append(_ref.euler)
        self.euler_ds.append(_ref.euler_d)
        self.euler_dds.append(_ref.euler_dd)
        self.rvel_bodys.append(_ref.rvel_body)
        self.raccel_bodys.append(_ref.raccel_body)


    def plot(self, time):
        eulers = np.array(self.eulers)
        euler_ds = np.array(self.euler_ds)
        euler_dds = np.array(self.euler_dds)
        plots = [("$\phi$",    "deg", np.rad2deg(eulers[:,0])),
                 ("$\\theta$", "deg", np.rad2deg(eulers[:,1])),
                 ("$\\psi$",   "deg", np.rad2deg(eulers[:,2])),
                 ("$\dot{\phi}$",    "deg/s", np.rad2deg(euler_ds[:,0])),
                 ("$\dot{\\theta}$", "deg/s", np.rad2deg(euler_ds[:,1])),
                 ("$\dot{\\psi}$",   "deg/s", np.rad2deg(euler_ds[:,2])),
                 ("$\ddot{\phi}$",    "deg/s2", np.rad2deg(euler_dds[:,0])),
                 ("$\ddot{\\theta}$", "deg/s2", np.rad2deg(euler_dds[:,1])),
                 ("$\ddot{\\psi}$",   "deg/s2", np.rad2deg(euler_dds[:,2]))]
        for i, (title, ylab, data) in enumerate(plots):
            ax = plt.subplot(3, 3, i+1)
            plt.plot(time, data)
            p3_pu.decorate(ax, title=title, ylab=ylab)

        plt.figure()
        rvel_bodys = np.array(self.rvel_bodys)
        raccel_bodys = np.array(self.raccel_bodys)
        plots = [("$\phi$",    "deg", np.rad2deg(eulers[:,0])),
                 ("$\\theta$", "deg", np.rad2deg(eulers[:,1])),
                 ("$\\psi$",   "deg", np.rad2deg(eulers[:,2])),
                 ("$p$",    "deg/s", np.rad2deg(rvel_bodys[:,0])),
                 ("$q$", "deg/s", np.rad2deg(rvel_bodys[:,1])),
                 ("$r$",   "deg/s", np.rad2deg(rvel_bodys[:,2])),
                 ("$\dot{p}$",    "deg/s2", np.rad2deg(raccel_bodys[:,0])),
                 ("$\dot{q}$",   "deg/s2", np.rad2deg(raccel_bodys[:,1])),
                 ("$\dot{r}$",   "deg/s2", np.rad2deg(raccel_bodys[:,2]))]
        for i, (title, ylab, data) in enumerate(plots):
            ax = plt.subplot(3, 3, i+1)
            plt.plot(time, data)
            p3_pu.decorate(ax, title=title, ylab=ylab)
        # test num integ
        dt = time[1]-time[0]
        ps = np.cumsum(raccel_bodys[:,0])*dt+rvel_bodys[0,0]
        plt.subplot(3,3,4); plt.plot(time, np.rad2deg(ps))
        qs = np.cumsum(raccel_bodys[:,1])*dt+rvel_bodys[0,1]
        plt.subplot(3,3,5); plt.plot(time, np.rad2deg(qs))
        rs = np.cumsum(raccel_bodys[:,2])*dt+rvel_bodys[0,2]
        plt.subplot(3,3,6); plt.plot(time, np.rad2deg(rs))
        

class PitchCtl:
    def __init__(self, dm, Xe, Ue, dt):
        self.Xe, self.Ue, self.dt = np.asarray(Xe), Ue, dt
        self.theta_e, self.q_e = Xe[p3_fr.SixDOFAeroEuler.sv_theta], Xe[p3_fr.SixDOFAeroEuler.sv_q] 
        self.sum_theta_err = 0
        self.ref = p3_u.SecOrdLinRef(omega=5., xi=0.9, sats=[np.deg2rad(45.), np.deg2rad(100.)])  # vel, accel
        self.compute_gain(dm, Xe, Ue, dt)
        
    def compute_gain(self, dm, Xe, Ue, dt):
        self.h_q, self.h_qd = -0.1, -0.01
        if 0: # empirical
            self.htkp1, self.hqkp1, self.htk, self.hqk = 0, 0, 0, 0
            #self.h_theta, self.h_q, self.h_qd = 0., -5.1, -2.
            self.k_theta, self.k_q, self.k_itheta = -6., -3., -0.075
        else: # learned
            self.htkp1, self.hqkp1, self.htk, self.hqk = -0.4237, -0.9071, -0.4234, 0.7884
            #self.k_theta, self.k_q, self.k_itheta  = -0.848, -0.053, -0.005
            self.k_theta, self.k_q, self.k_itheta  = -0.848, -0.053, -0.005
        print('theta loop gains')
        print('  feedback p:{} d:{} (i:{})'.format(self.k_theta, self.k_q, self.k_itheta))
        print('  feedforward {} {} {} {}'.format(self.htkp1, self.hqkp1, self.htk, self.hqk))
 

    def reset(self, theta, q):
        self.ref.reset(np.array([theta-self.theta_e, q-self.q_e, 0]))
        self.sum_theta_err = 0
     
    def get(self, X, t, theta_sp):
        sv = p3_fr.SixDOFAeroEuler
        theta_refk, q_refk = self.ref.X[:2]
        theta_refkp1, q_refkp1, qd_refkp1 = self.ref.run(self.dt, theta_sp-self.theta_e)
        # feedforward
        if 1:
            d_ele = self.htkp1*theta_refkp1 + self.hqkp1*q_refkp1 + self.htk*theta_refk + self.hqk*q_refk
        else:
            d_ele = 0.
        d_theta, d_q = X[sv.sv_theta]-self.theta_e, X[sv.sv_q]-self.q_e
        theta_err, q_err = d_theta-theta_refkp1, d_q-q_refkp1
        #theta_e, q_e = np.clip([theta_e, q_e], np.deg2rad([-30, -80]), np.deg2rad([30, 80]))
        self.sum_theta_err += theta_err

        # feedback
        d_ele -= self.k_theta*theta_err + self.k_q*q_err  + self.k_itheta*self.sum_theta_err
        return d_ele


    def get2(self, X, t, theta_sp, euler_ref):
        _s = p3_fr.SixDOFAeroEuler
        theta_err = X[_s.sv_theta] - euler_ref.theta
        q_err = X[_s.sv_q] - euler_ref.q
        self.sum_theta_err += theta_err
        # feedback
        d_ele = -(self.k_theta*theta_err + self.k_q*q_err  + self.k_itheta*self.sum_theta_err)
        # feedforward
        d_ele += self.h_q*euler_ref.rvel_body[1] + self.h_qd*euler_ref.raccel_body[1]
        return d_ele
    

class RollCtl:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = np.asarray(Xe), Ue
        self.dt = dt
        self.k_phi, self.k_p = -3.5, -0.75
        self.ref = p3_u.SecOrdLinRef(omega=6, xi=0.9, sats=[6., 50.])  # vel, accel
        self.h_p, self.h_pd = -0.15, -0.01
        
    def reset(self, phi, phi_d=0):
        self.ref.reset(np.array([phi, phi_d, 0]))
        
    def get(self, X, t, phi_sp, euler_ref):
        phi_ref, p_ref, pd_ref = self.ref.run(self.dt, phi_sp)
        phi_e, p_e = X[p3_fr.SixDOFAeroEuler.sv_phi]-euler_ref.phi, X[p3_fr.SixDOFAeroEuler.sv_p]-euler_ref.p
        #phi_e, p_e = X[p3_fr.SixDOFAeroEuler.sv_phi]-phi_ref, X[p3_fr.SixDOFAeroEuler.sv_p]-p_ref
        phi_e, p_e = np.clip([phi_e, p_e], np.deg2rad([-30, -150]), np.deg2rad([30, 150]))
        d_aile = -self.k_phi*phi_e -self.k_p*p_e
        # feedforward
        d_aile +=  self.h_p*euler_ref.rvel_body[0] + self.h_pd*euler_ref.raccel_body[0]
        return d_aile

class AttCtl:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = Xe, Ue
        self.dm = dm
        self.dt = dt
        self.euler_ref = EulerRef()
        self.pitch_ctl = PitchCtl(dm, Xe, Ue, dt)
        self.roll_ctl = RollCtl(Xe, Ue, dm, dt)
        self.sat_ail = np.array([np.deg2rad(-30), np.deg2rad(30)])
        self.sat_ele = np.array([np.deg2rad(-20), np.deg2rad(20)])
        self.k_beta, self.k_r = -0.2, -0.5
        self.h_r, self.h_rd = 0.1, 0.1
        
    def reset(self, t, X):
        sv = p3_fr.SixDOFAeroEuler
        phi, theta = X[sv.sv_phi], X[sv.sv_theta]
        p, q, r = X[sv.sv_slice_rvel]
        self.roll_ctl.reset(phi, p)
        self.pitch_ctl.reset(theta, q)
        self.euler_ref.reset(phi, theta, p, q)
        
    def get(self, t, X, phi_sp, theta_sp):
        va = X[p3_fr.SixDOFAeroEuler.sv_va]
        self.euler_ref.run(self.dt, phi_sp, theta_sp, va)
        U = np.array(self.Ue)
        U[self.dm.iv_da()] += self.roll_ctl.get( X, t, phi_sp, self.euler_ref)
        #U[self.dm.iv_de()] += self.pitch_ctl.get( X, t, theta_sp)
        U[self.dm.iv_de()] += self.pitch_ctl.get2( X, t, theta_sp, self.euler_ref)
        U[self.dm.iv_da()] = np.clip(U[self.dm.iv_da()], *self.sat_ail)
        U[self.dm.iv_de()] = np.clip(U[self.dm.iv_de()], *self.sat_ele)
        # yaw damper??
        beta, r_err = X[p3_fr.SixDOFAeroEuler.sv_beta], X[p3_fr.SixDOFAeroEuler.sv_r] - self.euler_ref.r
        U[self.dm.iv_dr()] += self.k_beta*beta + self.k_r*r_err
        U[self.dm.iv_dr()] += self.h_r * self.euler_ref.rvel_body[2] + self.h_rd * self.euler_ref.raccel_body[2]
        return U




#
# debug
#

    
class AttCtl2:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = Xe, Ue
        self.dm = dm

        
    def get(self, t, X, phi_sp, theta_sp):
        U = np.array(self.Ue)
        theta_err = X[p3_fr.SixDOFAeroEuler.sv_theta] - theta_sp
        q_err = X[p3_fr.SixDOFAeroEuler.sv_q]
        de = 1.*theta_err +0.2*q_err
        U[self.dm.iv_de()] += de
        return U



    # testing pitch control
       # if 0:
       #      A, B = dm.get_jacobian(Xe, Ue)
       #      sv = p3_fr.SixDOFAeroEuler; iv_ele = 2
       #      A1 = np.array([[A[sv.sv_alpha, sv.sv_alpha], A[sv.sv_alpha, sv.sv_theta], A[sv.sv_alpha, sv.sv_q]],
       #                     [A[sv.sv_theta, sv.sv_alpha], A[sv.sv_theta, sv.sv_theta], A[sv.sv_theta, sv.sv_q]],
       #                     [A[sv.sv_q,     sv.sv_alpha], A[sv.sv_q,     sv.sv_theta], A[sv.sv_q,     sv.sv_q]]])
       #      B1 = np.array([[B[sv.sv_alpha, iv_ele]],
       #                     [B[sv.sv_theta, iv_ele]],
       #                     [B[sv.sv_q,     iv_ele]]])
       #  if 0:  # pole placement
       #      eigva, eigve = np.linalg.eig(A1)
       #      om, xi = 9, 0.99; l1, l2 = p3_u.omxi_to_lambda(om, xi)
       #      print(' control objective cl poles: om {} xi {} l1 {} l2 {}'.format(om, xi, l1, l2))
       #      K1 = control.matlab.place(A1, B1, np.array([l1, l2, -om]))
       #      print('gain {}'.format(K1))
       #      # closed loop poles
       #      A1cl = A1 - np.dot(B1, K1)
       #      print('closed loop poles {}'.format(np.linalg.eig(A1cl)[0]))
       #      self.k_alpha, self.k_theta, self.k_q = K1.squeeze()
            
       #  if 0: # LQR
       #      Q, R = np.diag([1., 1., 0.1]), np.diag([1])
       #      (K2, X2, E) = control.matlab.lqr(A1, B1, Q, R)
       #      print('gain {}'.format(K2))
       #      # closed loop poles
       #      A1cl = A1 - np.dot(B1, K2)
       #      print('closed loop poles {}'.format(np.linalg.eig(A1cl)[0]))
       #      self.k_alpha, self.k_theta, self.k_q = K2.squeeze()
       #  #pdb.set_trace()
        
