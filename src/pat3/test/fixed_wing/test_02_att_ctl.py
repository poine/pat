#! /usr/bin/env python
import os, sys, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb


'''
  Testing Piloting loops
'''

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil
import pat3.utils as p3_u
import pat3.atmosphere as p3_atm
import pat3.frames as p3_fr
import pat3.plot_utils as p3_pu

def get_sim_defaults(t0=0, tf=5., dt=0.005, trim_args = {'h':0, 'va':12, 'gamma':0.}):
    trim_args = {'h':0, 'va':12, 'throttle':0., 'phi':np.deg2rad(0.)}
    dm = p1_fw_dyn.DynamicModel(param_filename)
    #dm = p1_fw_dyn.DynamicModel_ee(param_filename)
    #dm = p1_fw_dyn.DynamicModel_eq(param_filename) # no yet
    time = np.arange(t0, tf, dt)
    Xe, Ue = dm.trim(trim_args, report=True, debug=False)
    Xe_ae = dm.state_as_six_dof_euclidian_euler(Xe)
    phi_sp = np.ones(len(time))*Xe_ae[p3_fr.SixDOFEuclidianEuler.sv_phi]
    theta_sp = np.ones(len(time))*Xe_ae[p3_fr.SixDOFEuclidianEuler.sv_theta]

    return dm, time, Xe, Ue, phi_sp, theta_sp

def run_simulation(dm, time, Xe, Ue, X0, phi_sp, theta_sp, plot=True):
    atm = p3_atm.AtmosphereCstWind([0, 0, 0])
    X = np.zeros((len(time), dm.sv_size))
    X_act = np.zeros((len(time), dm.input_nb()))
    theta_ref1 = np.zeros((len(time), 3))
    phi_ref =  np.zeros(len(time))
    theta_ref =  np.zeros(len(time))
    p_ref, q_ref, r_ref = np.zeros(len(time)), np.zeros(len(time)), np.zeros(len(time)) 
    U = np.zeros((len(time),  dm.input_nb()))
    X_act[0] = Ue
    X[0] = dm.reset(X0, time[0], X_act0=Ue)
    ctl = p3_pil.AttCtl(Xe, Ue, dm, time[1]-time[0])
    ctl.reset(time[0], X[0])
    #pdb.set_trace()
    theta_ref1[0] = ctl.pitch_ctl.ref.X
    #ctl = p3_pil.AttCtl2(Xe, Ue, dm, time[1]-time[0])
    for i in range(1, len(time)):
        U[i-1] = ctl.get(time[i-1], X[i-1], phi_sp[i-1], theta_sp[i-1])
        theta_ref1[i-1] = ctl.pitch_ctl.ref.X
        phi_ref[i-1] = ctl.euler_ref.phi
        theta_ref[i-1] = ctl.euler_ref.theta
        p_ref[i-1], q_ref[i-1], r_ref[i-1] = ctl.euler_ref.p, ctl.euler_ref.q, ctl.euler_ref.r
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1], atm=atm)
        X_act[i] = dm.X_act

    U[-1] = ctl.get(time[-1], X[-1], phi_sp[-1], theta_sp[-1])
    phi_ref[-1] = ctl.euler_ref.phi
    theta_ref[-1] = ctl.euler_ref.theta
    if plot:
        dm.plot_trajectory(time, X, X_act)

        plt.subplot(5,3,7) # phi
        plt.plot(time, np.rad2deg(phi_sp), label='sp')
        plt.plot(time, np.rad2deg(phi_ref), label='ref')
        plt.legend()

        plt.subplot(5,3,8) # theta
        plt.plot(time, np.rad2deg(theta_sp), label='setpoint')
        #plt.plot(time, np.rad2deg(theta_ref1[:,0]+ctl.pitch_ctl.theta_e), label='reference')
        plt.plot(time, np.rad2deg(theta_ref), label='reference2')
        plt.legend()

        plt.subplot(5,3,10) # p
        plt.plot(time, np.rad2deg(p_ref), label='reference2')
        plt.legend()

        plt.subplot(5,3,11) # q
        #plt.plot(time, np.rad2deg(theta_ref1[:,1]+ctl.pitch_ctl.q_e), label='reference')
        plt.plot(time, np.rad2deg(q_ref), label='reference2')
        plt.legend()

        plt.subplot(5,3,12) # r
        plt.plot(time, np.rad2deg(r_ref), label='reference2')
        plt.legend()
        
        plt.subplot(5,3,13); plt.plot(time, 100*U[:,dm.iv_dth()], label='input')     # throttle
        plt.subplot(5,3,14); plt.plot(time, np.rad2deg(U[:,dm.iv_da()]), alpha=0.5)  # aileron
        plt.subplot(5,3,15); plt.plot(time, np.rad2deg(U[:,dm.iv_de()]), alpha=0.5)  # elevator
    return time, X, U

def test_pert_theta(ampl=np.deg2rad(-1.)):
    dm, time, Xe, Ue, phi_sp, theta_sp = get_sim_defaults(tf=10.)
    X0 = np.array(Xe); X0[p3_fr.SixDOFEuclidianEuler.sv_theta] += ampl
    return run_simulation(dm, time, Xe, Ue, X0, phi_sp, theta_sp, plot=True)


def test_step_phi(ampl=np.deg2rad(20.), _dur=10., _p=5.):
    dm, time, Xe, Ue, phi_sp, theta_sp = get_sim_defaults(tf=_dur)
    phi_sp = np.array([p3_u.step(_t, a=ampl, p=_p, dt=0) for _t in time])
    return run_simulation(dm, time, Xe, Ue, Xe, phi_sp, theta_sp, plot=True)

def test_step_theta(ampl=np.deg2rad(1.), _dur=10., _p=10.):
    dm, time, Xe, Ue, phi_sp, theta_sp = get_sim_defaults(tf=_dur)
    theta_sp = Xe[p3_fr.SixDOFAeroEuler.sv_theta]+np.array([p3_u.step(_t, a=ampl, p=_p, dt=_p/4.) for _t in time])
    return run_simulation(dm, time, Xe, Ue, Xe, phi_sp, theta_sp, plot=True)

def test_step_phi_theta(phi_ampl=np.deg2rad(20.), theta_ampl=np.deg2rad(1), _dur=10., _p=5.):
    dm, time, Xe, Ue, phi_sp, theta_sp = get_sim_defaults(tf=_dur)
    phi_sp = np.array([p3_u.step(_t, a=phi_ampl, p=_p, dt=0) for _t in time])
    theta_sp = Xe[p3_fr.SixDOFAeroEuler.sv_theta]+np.array([p3_u.step(_t, a=theta_ampl, p=10., dt=2.5) for _t in time])
    return run_simulation(dm, time, Xe, Ue, Xe, phi_sp, theta_sp, plot=True)
    
    
def main(param_filename):
    #test_pert_theta()
    #test_step_phi(ampl=np.deg2rad(45.), _dur=20., _p=10)
    test_step_theta(ampl=np.deg2rad(5.), _dur=90., _p=70. )
    #test_step_phi_theta(phi_ampl=np.deg2rad(20.), theta_ampl=np.deg2rad(5.), _dur=10., _p=5.)
    plt.show()  
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = p3_u.pat_ressource('data/vehicles/cularis.xml')
    #param_filename = p3_u.pat_ressource('data/vehicles/funjet_avl.xml') # not yet :(
    main(param_filename)
