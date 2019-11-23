#! /usr/bin/env python
import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil
import pat3.vehicles.fixed_wing.guidance as p3_guid

import pat3.atmosphere as p3_atm
import pat3.frames as p3_fr
import pat3.trajectory_3D as p3_traj3d
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.algebra as p3_alg

import control.matlab

import test_02_att_ctl


#def norm_mpi_pi(v): return ( v + np.pi) % (2 * np.pi ) - np.pi



    
def run_simulation(dm, ref_traj, tf=60.5, dt=0.01, trim_args = {'h':0, 'va':12, 'gamma':0}, plot=False, atm=None):
    time = np.arange(0, tf, dt)
    X = np.zeros((len(time), dm.sv_size))
    U = np.zeros((len(time),  dm.input_nb()))
    carrots, att_sp = np.zeros((len(time),  3)), np.zeros((len(time), 2))
    ctl = p3_guid.Guidance(dm, ref_traj, trim_args, dt)
    X[0] = dm.reset(ctl.Xe, t0=time[0], X_act0=None)#ctl.Ue)
    for i in range(1, len(time)):
        Xee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(X[i-1])
        U[i-1] = ctl.get(time[i-1], X[i-1], Xee)
        carrots[i-1] = ctl.carrot; att_sp[i-1] = ctl.phi_sp, ctl.theta_sp 
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1], atm)
    U[-1] = ctl.get(time[-1], X[-1]); carrots[-1] = ctl.carrot; att_sp[-1] = ctl.phi_sp, ctl.theta_sp 
    if plot:
        dm.plot_trajectory(time, X, U)
        plt.subplot(5,3,1); plt.plot(time, carrots[:,0])
        plt.subplot(5,3,2); plt.plot(time, carrots[:,1])
        plt.subplot(5,3,3); plt.plot(time, carrots[:,2])
        plt.subplot(5,3,7); plt.plot(time, np.rad2deg(att_sp[:,0]), label='setpoint')
        plt.subplot(5,3,8); plt.plot(time, np.rad2deg(att_sp[:,1]), label='setpoint')
        plt.show()  
    return time, X, U


    
    
def main(param_filename, trim_args = {'h':0, 'va':11, 'gamma':0}):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    ref_traj = p3_traj3d.CircleRefTraj(c=[0, 0, 0], r=20)
    #time, X, U = test_02_att_ctl.run_simulation(dm, plot=False)
    time, X, U = run_simulation(dm, ref_traj, plot=True)
    
    #atm =  p3_atm.Atmosphere()
    p3_pu.plot_3D_traj(ref_traj, X)
    #p3_pu.plot_3D_wind(atm)
    plt.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    dirname, filename = os.path.split(os.path.abspath(__file__))
    param_filename = os.path.abspath(os.path.join(dirname, '../../../../data/vehicles/cularis.xml'))
    main(param_filename)
