#! /usr/bin/env python
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

''' Testing fixed wing in open loop '''

import pat3.utils as p3_u
import pat3.frames as p3_fr
import pat3.dynamics as pat_dyn
import pat3.vehicles.fixed_wing.simple_6dof_fdm as fw_dyn
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn

DEFAULT_AC='cularis'
#DEFAULT_AC='skywalker_x8'
def get_default_dm(ac_name=DEFAULT_AC):
    #return p1_fw_dyn.DynamicModel(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))
    return p1_fw_dyn.DynamicModel_ee(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))

def get_sim_defaults(t0=0, tf=5., dt=0.005, trim_args = {'h':0, 'va':10, 'gamma':0}, ac_name='cularis'):
    dm = get_default_dm(ac_name)
    time = np.arange(t0, tf, dt)
    Xe, Ue = dm.trim(trim_args, report=True)
    return dm, time, Xe, Ue

def run_simulation(dm, time, X0, X_act0, U, atm=None):
    X = np.zeros((len(time), dm.sv_size))
    X_act = np.zeros((len(time), dm.input_nb()))
    X[0] = dm.reset(X0, time[0], X_act0)
    X_act[0] = dm.X_act
    for i in range(1, len(time)):
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1], atm=atm)
        X_act[i] = dm.X_act
    return time, X, X_act

def sim_pert(pert_axis, pert_val, label, t0=0, tf=1., dt=0.01, trim_args = {'h':0, 'va':12, 'gamma':0}, plot=True):
    dm, time, Xe, Ue = get_sim_defaults(t0, tf, dt, trim_args)
    dX0 = np.zeros(dm.sv_size)
    U = Ue*np.ones((len(time), dm.input_nb()))
    Xpert = np.zeros(dm.sv_size); Xpert[pert_axis] = pert_val
    X0, X_act0 = Xe+Xpert, Ue # we start with equilibred actuators
    time, X, X_act = run_simulation(dm, time, X0, X_act0, U, atm=None)
    if plot:
        dm.plot_trajectory(time, X, X_act, window_title=label)
        plt.subplot(5,3,7); plt.plot(time, np.rad2deg(Xe[dm.sv_phi]*np.ones(len(time))))
        plt.subplot(5,3,8); plt.plot(time, np.rad2deg(Xe[dm.sv_theta]*np.ones(len(time))))
    return time, X, U

def sim_pert_p(): return sim_pert(p3_fr.SixDOFAeroEuler.sv_p, np.deg2rad(10.), 'perturbation roll rate', tf=10.)
def sim_pert_q(): return sim_pert(p3_fr.SixDOFAeroEuler.sv_q, np.deg2rad(10.), 'perturbation pitch rate', tf=5.)
def sim_pert_r(): return sim_pert(p3_fr.SixDOFAeroEuler.sv_r, np.deg2rad(10.), 'perturbation yaw rate', tf=10.)
def sim_pert_phi(): return sim_pert(p3_fr.SixDOFAeroEuler.sv_phi, np.deg2rad(10.), 'perturbation roll', tf=10.)
def sim_pert_theta(): return sim_pert(p3_fr.SixDOFAeroEuler.sv_theta, np.deg2rad(2.), 'perturbation pitch', tf=45.)
def sim_pert_va(): return sim_pert(p3_fr.SixDOFAeroEuler.sv_va, 1., 'perturbation air vel', tf=40.)             # 40s to see phugoid
def sim_pert_alpha(): return sim_pert(p3_fr.SixDOFAeroEuler.sv_alpha, np.deg2rad(1.), 'perturbation alpha', tf=10.) # 
def sim_pert_beta(): return sim_pert(p3_fr.SixDOFAeroEuler.sv_beta, np.deg2rad(1.), 'perturbation beta', tf=10.)  # 


def sim_step_act(act_id, act_val, label, t0=0., tf=10., dt=0.01, trim_args = {'h':0, 'va':12, 'gamma':0}, plot=True):
    dm, time, Xe, Ue = get_sim_defaults(t0, tf, dt, trim_args)
    U = Ue*np.ones((len(time), dm.input_nb()))
    U[:,act_id] += p3_u.step_vec(time, a=act_val, p=10., dt=2.5)
    X0, X_act0 = Xe, Ue
    time, X, X_act = run_simulation(dm, time, Xe, Ue, U, atm=None)
    if plot:
        dm.plot_trajectory(time, X, X_act, window_title=label, label=dm.get_name()) 
    return time, X, U

def sim_step_thr(): return sim_step_act(0, 0.05, 'step throttle', tf=20.)
def sim_step_ail(): return sim_step_act(1, np.deg2rad(1.), 'step aileron')
def sim_step_ele(): return sim_step_act(2, np.deg2rad(0.5), 'step elevator', tf=20.)
def sim_step_rud(): return sim_step_act(3, np.deg2rad(10.),   'step rudder')
def sim_step_flap(): return sim_step_act(4, np.deg2rad(0.5), 'step flap')

                
def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)

    if 0:
        sim_pert_p()
        sim_pert_q()
        sim_pert_r()
    if 0:
        sim_pert_phi()
        sim_pert_theta()
        sim_pert_va()
        sim_pert_alpha()
        sim_pert_beta()
    if 1:
        #sim_step_thr()
        sim_step_ele()
        #sim_step_ail()
        #sim_step_rud()
        #sim_step_flap()
    plt.show()


    
if __name__ == "__main__":
    main()
