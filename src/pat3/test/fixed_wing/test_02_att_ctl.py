#! /usr/bin/env python
import os, sys, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import control.matlab

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil
import pat3.utils as p3_u
import pat3.atmosphere as p3_atm
import pat3.frames as p3_fr
import pat3.plot_utils as p3_pu

def get_sim_defaults(dm, t0=0, tf=5., dt=0.005, trim_args = {'h':0, 'va':12, 'gamma':0}):
    time = np.arange(t0, tf, dt)
    Xe, Ue = dm.trim(trim_args, report=True, debug=False)
    phi_sp = np.ones(len(time))*Xe[dm.sv_phi]
    theta_sp = np.ones(len(time))*Xe[dm.sv_theta]
    return time, Xe, Ue, phi_sp, theta_sp

def run_simulation(dm, time, Xe, Ue, phi_sp, theta_sp, plot=True):
    atm = p3_atm.AtmosphereCstWind([0, 0, 0])
    X = np.zeros((len(time), dm.sv_size()))
    X_act = np.zeros((len(time), dm.input_nb()))
    U = np.zeros((len(time),  dm.input_nb()))
    X[0] = dm.reset(Xe)
    ctl = p3_pil.AttCtl(Xe, Ue, dm, time[1]-time[0])
    for i in range(1, len(time)):
        U[i-1] = ctl.get(time[i-1], X[i-1], phi_sp[i-1], theta_sp[i-1])
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1], atm=atm)
        X_act[i] = dm.X_act
    U[-1] = ctl.get(time[-1], X[-1], phi_sp[-1], theta_sp[-1])
    if plot:
        dm.plot_trajectory(time, X, X_act)
        plt.subplot(5,3,7); plt.plot(time, np.rad2deg(phi_sp))
        plt.subplot(5,3,8); plt.plot(time, np.rad2deg(theta_sp))
        plt.subplot(5,3,13); plt.plot(time, 100*U[:,dm.iv_dth()]) # throttle
        plt.subplot(5,3,14); plt.plot(time, np.rad2deg(U[:,dm.iv_da()]), alpha=0.5)  # aileron
        plt.subplot(5,3,15); plt.plot(time, np.rad2deg(U[:,dm.iv_de()]), alpha=0.5)  # elevator
    return time, X, U


    

def test_step_phi(dm, ampl=np.deg2rad(20.)):
    time, Xe, Ue, phi_sp, theta_sp = get_sim_defaults(dm)
    phi_sp = np.array([p3_u.step(_t, a=ampl, p=5., dt=1.25) for _t in time])
    return run_simulation(dm, time, Xe, Ue, phi_sp, theta_sp, plot=True)

def test_step_theta(dm, ampl=np.deg2rad(1.)):
    time, Xe, Ue, phi_sp, theta_sp = get_sim_defaults(dm)
    theta_sp = Xe[p3_fr.SixDOFAeroEuler.sv_theta]+np.array([p3_u.step(_t, a=ampl, p=10., dt=2.5) for _t in time])
    return run_simulation(dm, time, Xe, Ue, phi_sp, theta_sp, plot=True)
    
    
def main(param_filename):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    #test_step_phi(dm)
    test_step_theta(dm)
    plt.show()  
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
    main(param_filename)
