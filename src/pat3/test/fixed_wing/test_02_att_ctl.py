#! /usr/bin/env python
import os, sys, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import control.matlab

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.utils as p3_u

class PitchCtl:
    def __init__(self, Xe, Ue, dm, dt):
        self.Xe, self.Ue = np.asarray(Xe), Ue
        self.dt = dt
        A, B = dm.get_jacobian(Xe, Ue)
        # if 0:
        #     A1, B1 = A[dm.sv_phi:, dm.sv_phi:], B[dm.sv_phi:, 2]
        #     pdb.set_trace()
        #     control.matlab.place(A1, B1, [])
        self.k_itheta, self.k_theta, self.k_q = -0.075, -20., -1.5
        self.sum_theta_err = 0
        self.h_theta, self.h_q, self.h_qd = -0.1, -0.1, -0.
        self.ref = p3_u.SecOrdLinRef(omega=6, xi=0.9, sats=[6., 50.])  # vel, accel
        
    def get(self, X, t, theta_sp):
        theta_ref, q_ref, qd_ref = self.ref.run(self.dt, theta_sp)
        theta_e, q_e = X[p1_fw_dyn.sv_theta]-theta_ref, X[p1_fw_dyn.sv_q]-q_ref
        #theta_e, q_e = np.clip([theta_e, q_e], np.deg2rad([-30, -80]), np.deg2rad([30, 80]))
        self.sum_theta_err += theta_e
        # feedback
        d_ele = -self.k_theta*theta_e -self.k_q*q_e -self.k_itheta*self.sum_theta_err
        # feedforward
        d_ele += self.h_theta*theta_ref + self.h_q*q_ref + self.h_qd*qd_ref 
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
    
    
def run_simulation(dm, tf=40.5, dt=0.01, trim_args = {'h':0, 'va':12, 'gamma':0}, plot=True):
    Xe, Ue = dm.trim(trim_args, debug=True)
    time = np.arange(0, tf, dt)
    X = np.zeros((len(time), dm.sv_size))
    X_act = np.zeros((len(time), dm.input_nb()))
    U = np.zeros((len(time),  dm.input_nb()))
    #phi_sp = np.array([p3_u.step(_t, a=np.deg2rad(1.), p=10., dt=5.) for _t in time])
    phi_sp = np.array([p3_u.step(_t, a=np.deg2rad(20.), p=5., dt=1.25) for _t in time])
    #phi_sp = np.zeros(len(time))
    theta_sp = np.array([p3_u.step(_t, a=np.deg2rad(1.), p=40., dt=0.) for _t in time]) + np.ones(len(time))*Xe[dm.sv_theta]
    #theta_sp = np.ones(len(time))*Xe[dm.sv_theta]
    X[0] = dm.reset(Xe)
    ctl = AttCtl(Xe, Ue, dm, dt)
    for i in range(1, len(time)):
        U[i-1] = ctl.get(time[i-1], X[i-1], phi_sp[i-1], theta_sp[i-1])
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1])
        X_act[i] = dm.X_act
    U[-1] = ctl.get(time[-1], X[-1], phi_sp[-1], theta_sp[-1])
    if plot:
        dm.plot_trajectory(time, X, X_act)
        plt.subplot(5,3,7); plt.plot(time, np.rad2deg(phi_sp))
        plt.subplot(5,3,8); plt.plot(time, np.rad2deg(theta_sp))
        plt.subplot(5,3,13); plt.plot(time, 100*U[:,dm.iv_dth()]) # throttle
        plt.subplot(5,3,14); plt.plot(time, np.rad2deg(U[:,dm.iv_da()]), alpha=0.2)  # aileron
        plt.subplot(5,3,15); plt.plot(time, np.rad2deg(U[:,dm.iv_de()]), alpha=0.2)  # elevator
        plt.show()  
    return time, X, U
        
def main(param_filename):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    run_simulation(dm)
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    dirname, filename = os.path.split(os.path.abspath(__file__))
    param_filename = os.path.abspath(os.path.join(dirname, '../../../../data/vehicles/cularis.xml'))
    main(param_filename)
