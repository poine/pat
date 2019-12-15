import os, sys, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import control.matlab

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.utils as p3_u
import pat3.frames as p3_fr
import pat3.atmosphere as p3_atm

class PitchCtl:
    def __init__(self, dm, Xe, Ue, dt):
        self.Xe, self.Ue, self.dt = np.asarray(Xe), Ue, dt
        #self.k_itheta, self.k_theta, self.k_q = -0.075, -10., -1.5
        self.sum_theta_err = 0
        self.ref = p3_u.SecOrdLinRef(omega=5, xi=0.9, sats=[np.deg2rad(45.), np.deg2rad(100.)])  # vel, accel
        self.compute_gain(dm, Xe, Ue, dt)
        
    def compute_gain(self, dm, Xe, Ue, dt):
        if 1:  # empirical
            #self.h_theta, self.h_q, self.h_qd = -0.001, -0.001, 0.
            self.h_theta, self.h_q, self.h_qd = 0., -5.1, -2.
            #self.h_theta, self.h_q, self.h_qd = 0., 0., 0.
            self.k_alpha, self.k_itheta, self.k_theta, self.k_q = 0., -0.075, -6., -3.
            #self.k_alpha, self.k_itheta, self.k_theta, self.k_q = 0., 0., -1., 0.
            print('theta loop gains {} {} {} {}'.format(self.k_alpha, self.k_theta, self.k_q, self.k_itheta))
        if 0:
            A, B = dm.get_jacobian(Xe, Ue)
            sv = p3_fr.SixDOFAeroEuler; iv_ele = 2
            A1 = np.array([[A[sv.sv_alpha, sv.sv_alpha], A[sv.sv_alpha, sv.sv_theta], A[sv.sv_alpha, sv.sv_q]],
                           [A[sv.sv_theta, sv.sv_alpha], A[sv.sv_theta, sv.sv_theta], A[sv.sv_theta, sv.sv_q]],
                           [A[sv.sv_q,     sv.sv_alpha], A[sv.sv_q,     sv.sv_theta], A[sv.sv_q,     sv.sv_q]]])
            B1 = np.array([[B[sv.sv_alpha, iv_ele]],
                           [B[sv.sv_theta, iv_ele]],
                           [B[sv.sv_q,     iv_ele]]])
        if 0:  # pole placement
            eigva, eigve = np.linalg.eig(A1)
            om, xi = 9, 0.99; l1, l2 = p3_u.omxi_to_lambda(om, xi)
            print(' control objective cl poles: om {} xi {} l1 {} l2 {}'.format(om, xi, l1, l2))
            K1 = control.matlab.place(A1, B1, np.array([l1, l2, -om]))
            print('gain {}'.format(K1))
            # closed loop poles
            A1cl = A1 - np.dot(B1, K1)
            print('closed loop poles {}'.format(np.linalg.eig(A1cl)[0]))
            self.k_alpha, self.k_theta, self.k_q = K1.squeeze()
            
        if 0: # LQR
            Q, R = np.diag([1., 1., 0.1]), np.diag([1])
            (K2, X2, E) = control.matlab.lqr(A1, B1, Q, R)
            print('gain {}'.format(K2))
            # closed loop poles
            A1cl = A1 - np.dot(B1, K2)
            print('closed loop poles {}'.format(np.linalg.eig(A1cl)[0]))
            self.k_alpha, self.k_theta, self.k_q = K2.squeeze()
        #pdb.set_trace()
        


        
    def get(self, X, t, theta_sp):
        sv = p3_fr.SixDOFAeroEuler
        alpha_ref = self.Xe[sv.sv_alpha]
        theta_ref, q_ref, qd_ref = self.ref.run(self.dt, theta_sp)
        alpha_e, theta_e, q_e = X[sv.sv_alpha]-alpha_ref, X[sv.sv_theta]-theta_ref, X[sv.sv_q]-q_ref
        #theta_e, q_e = np.clip([theta_e, q_e], np.deg2rad([-30, -80]), np.deg2rad([30, 80]))
        self.sum_theta_err += theta_e
        # feedforward
        d_ele = self.h_theta*theta_ref + self.h_q*q_ref + self.h_qd*qd_ref 
        # feedback
        d_ele -= self.k_alpha*alpha_e + self.k_theta*theta_e + self.k_q*q_e  + self.k_itheta*self.sum_theta_err
        return d_ele


class RollCtl:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = np.asarray(Xe), Ue
        self.dt = dt
        self.k_phi, self.k_p = -3.5, -0.75
        self.ref = p3_u.SecOrdLinRef(omega=6, xi=0.9, sats=[6., 50.])  # vel, accel
        
        
    def get(self, X, t, phi_sp):
        phi_ref, p_ref, pd_ref = self.ref.run(self.dt, phi_sp)
        phi_e, p_e = X[p3_fr.SixDOFAeroEuler.sv_phi]-phi_ref, X[p3_fr.SixDOFAeroEuler.sv_p]-p_ref
        phi_e, p_e = np.clip([phi_e, p_e], np.deg2rad([-30, -150]), np.deg2rad([30, 150]))
        d_aile = -self.k_phi*phi_e -self.k_p*p_e
        return d_aile

class AttCtl:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = Xe, Ue
        self.dm = dm
        self.pitch_ctl = PitchCtl(dm, Xe, Ue, dt)
        self.roll_ctl = RollCtl(Xe, Ue, dm, dt)
        self.sat_ail = np.array([np.deg2rad(-30), np.deg2rad(30)])
        self.sat_ele = np.array([np.deg2rad(-20), np.deg2rad(20)])
        
    def get(self, t, X, phi_sp, theta_sp):
        U = np.array(self.Ue)
        U[self.dm.iv_da()] += self.roll_ctl.get( X, t, phi_sp)
        U[self.dm.iv_de()] += self.pitch_ctl.get( X, t, theta_sp)
        U[self.dm.iv_da()] = np.clip(U[self.dm.iv_da()], *self.sat_ail)
        U[self.dm.iv_de()] = np.clip(U[self.dm.iv_de()], *self.sat_ele)
        return U
