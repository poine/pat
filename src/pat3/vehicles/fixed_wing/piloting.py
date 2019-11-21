import os, sys, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import control.matlab

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.utils as p3_u
import pat3.atmosphere as p3_atm

class PitchCtl:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = np.asarray(Xe), Ue
        self.dt = dt
        A, B = dm.get_jacobian(Xe, Ue)
        # if 0:
        #     A1, B1 = A[dm.sv_phi:, dm.sv_phi:], B[dm.sv_phi:, 2]
        #     pdb.set_trace()
        #     control.matlab.place(A1, B1, [])
        #self.k_itheta, self.k_theta, self.k_q = -0.075, -10., -1.5
        self.k_itheta, self.k_theta, self.k_q = -0.075, -6., -3.
        self.sum_theta_err = 0
        #self.h_theta, self.h_q, self.h_qd = -0.001, -0.001, 0.
        self.h_theta, self.h_q, self.h_qd = 0., -5.1, -1.
        self.ref = p3_u.SecOrdLinRef(omega=5, xi=0.9, sats=[np.deg2rad(45.), np.deg2rad(100.)])  # vel, accel
        
    def get(self, X, t, theta_sp):
        theta_ref, q_ref, qd_ref = self.ref.run(self.dt, theta_sp)
        theta_e, q_e = X[p1_fw_dyn.sv_theta]-theta_ref, X[p1_fw_dyn.sv_q]-q_ref
        #theta_e, q_e = np.clip([theta_e, q_e], np.deg2rad([-30, -80]), np.deg2rad([30, 80]))
        self.sum_theta_err += theta_e
        # feedforward
        d_ele = self.h_theta*theta_ref + self.h_q*q_ref + self.h_qd*qd_ref 
        # feedback
        d_ele -= self.k_theta*theta_e + self.k_q*q_e  + self.k_itheta*self.sum_theta_err
        return d_ele


class RollCtl:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = np.asarray(Xe), Ue
        self.dt = dt
        self.k_phi, self.k_p = -3.5, -0.75
        self.ref = p3_u.SecOrdLinRef(omega=6, xi=0.9, sats=[6., 50.])  # vel, accel
        
        
    def get(self, X, t, phi_sp):
        phi_ref, p_ref, pd_ref = self.ref.run(self.dt, phi_sp)
        phi_e, p_e = X[p1_fw_dyn.sv_phi]-phi_ref, X[p1_fw_dyn.sv_p]-p_ref
        phi_e, p_e = np.clip([phi_e, p_e], np.deg2rad([-30, -150]), np.deg2rad([30, 150]))
        d_aile = -self.k_phi*phi_e -self.k_p*p_e
        return d_aile

class AttCtl:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = Xe, Ue
        self.dm = dm
        self.pitch_ctl = PitchCtl(Xe, Ue, dm, dt)
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
