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

def get_sim_defaults(t0=0, tf=5., dt=0.005, trim_args = {'h':0, 'va':10, 'gamma':0}):
    dm = p1_fw_dyn.DynamicModel(os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml'))
    time = np.arange(t0, tf, dt)
    Xe, Ue = dm.trim(trim_args, debug=True)
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

def sim_pert_p(): sim_pert(p3_fr.SixDOFAeroEuler.sv_p, np.deg2rad(10.), 'perturbation roll rate', tf=10.)
def sim_pert_q(): sim_pert(p3_fr.SixDOFAeroEuler.sv_q, np.deg2rad(10.), 'perturbation pitch rate', tf=5.)
def sim_pert_r(): sim_pert(p3_fr.SixDOFAeroEuler.sv_r, np.deg2rad(10.), 'perturbation yaw rate', tf=10.)
def sim_pert_phi(): sim_pert(p3_fr.SixDOFAeroEuler.sv_phi, np.deg2rad(10.), 'perturbation roll', tf=10.)
def sim_pert_theta(): sim_pert(p3_fr.SixDOFAeroEuler.sv_theta, np.deg2rad(2.), 'perturbation pitch', tf=15.)
def sim_pert_va(): sim_pert(p3_fr.SixDOFAeroEuler.sv_va, 1., 'perturbation air vel', tf=40.)             # 40s to see phugoid
def sim_pert_alpha(): sim_pert(p3_fr.SixDOFAeroEuler.sv_alpha, np.deg2rad(1.), 'perturbation alpha', tf=10.) # 
def sim_pert_beta(): sim_pert(p3_fr.SixDOFAeroEuler.sv_beta, np.deg2rad(1.), 'perturbation beta', tf=10.)  # 


def sim_step_act(act_id, act_val, label, t0=0., tf=10., dt=0.01, trim_args = {'h':0, 'va':12, 'gamma':0}, plot=True):
    dm, time, Xe, Ue = get_sim_defaults(t0, tf, dt, trim_args)
    U = Ue*np.ones((len(time), dm.input_nb()))
    U[:,act_id] += p3_u.step_vec(time, a=act_val, p=10., dt=2.5)
    X0, X_act0 = Xe, Ue
    time, X, X_act = run_simulation(dm, time, Xe, Ue, U, atm=None)
    if plot:
        dm.plot_trajectory(time, X, X_act, window_title=label) 

def sim_step_thr(): sim_step_act(0, 0.01, 'step throttle', tf=20.)
def sim_step_ail(): sim_step_act(1, np.deg2rad(1.), 'step aileron')
def sim_step_ele(): sim_step_act(2, np.deg2rad(0.5), 'step elevator')
def sim_step_rud(): sim_step_act(3, np.deg2rad(0.5), 'step rudder')
def sim_step_flap(): sim_step_act(4, np.deg2rad(0.5), 'step flap')

def get_trim_defaults(trim_args={'h':0, 'va':12, 'gamma':0}):
    dm = p1_fw_dyn.DynamicModel(os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml'))
    Xe, Ue = dm.trim(trim_args, debug=True)
    return dm, Xe, Ue

def print_trim(trim_args={'h':0, 'va':12, 'gamma':0}):
    dm, Xe, Ue = get_trim_defaults(trim_args)
    print(dm.state_str(Xe))

def print_trim_dyn(trim_args={'h':0, 'va':12, 'gamma':0}):
    dm, Xe, Ue = get_trim_defaults(trim_args)
    A, B = dm.get_jacobian(Xe, Ue)
    _eval, _evel = np.linalg.eig(A)
    print('Modes:')
    for i, _ev in enumerate(_eval):
        if np.isreal(_ev):
            txt = 'lambda_{}  : {:.3f}'.format(i+1, np.real(_ev))
            if np.abs(np.real(_ev))>1e-6: txt += ' (tau {:.1f}s)'.format(np.abs(1./np.real(_ev)))
            print(txt)
        else:
            if i>0 and np.imag(_ev)==-np.imag(_eval[i-1]):
                txt = ''
            else:
                om = np.linalg.norm(_ev); xi = -np.real(_ev)/om
                txt = 'lambda_{}/{}: {:.3f} (om {:.2f} rad/s, xi {:.2f})'.format(i+1, i+2, _ev, om, xi)
                print(txt)
            
        
def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    #test_equilibrium() # for now broken fdm
    #print_trim()
    #print_trim_dyn()
    #sim_pert_p()
    #sim_pert_q()
    #sim_pert_r()
    #sim_pert_phi()
    #sim_pert_theta()
    #sim_pert_va()
    #sim_pert_alpha()
    #sim_pert_beta()
    #sim_step_thr()
    sim_step_ele()
    sim_step_ail()
    #sim_step_rud()
    #sim_step_flap()
    plt.show()


    
if __name__ == "__main__":
    main()
