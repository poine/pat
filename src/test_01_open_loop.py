#! /usr/bin/env python
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.algebra as pal

''' Testing quadrotor in open loop '''

def test_dyn(P):
    Xe, Ue = fdm.trim(P)
    dU = 0.001*np.array([1, -1, -1, 1])
    Xk, Uk = Xe, Ue
    phi0, theta0, psi0 = np.deg2rad(1.), 0., 0.
    Xk[fdm.sv_slice_quat] = pal.quat_of_euler([phi0, theta0, psi0])
    Xkp1 = fdm.disc_dyn(Xk, 0, Uk, 0.01, P)
    pdb.set_trace()

def step(t, a=-1., p=10., dt=0.): return a if math.fmod(t+dt, p) > p/2 else -a

def default_sim_env(P, tf=10, dt=0.01):
    time = np.arange(0, tf, dt) 
    Xe, Ue = fdm.trim(P)
    return time, Xe, Ue*np.ones((len(time), fdm.iv_size))
    
def sim_open_loop(P, time, X0, U):
    X = np.zeros((len(time), fdm.sv_size))
    X[0] = X0
    for i in range(1, len(time)):
        X[i] = fdm.disc_dyn(X[i-1], time[i-1], U[i-1], time[i]-time[i-1], P) 
    return X


def sim_pert_phi(P):
    time, Xe, U = default_sim_env(P)
    X0 = np.array(Xe)
    phi0, theta0, psi0 = np.deg2rad(1.), 0., 0.
    X0[fdm.sv_slice_quat] = pal.quat_of_euler([phi0, theta0, psi0])
    # pdb.set_trace()
    return time, sim_open_loop(P, time, X0, U), U, 'pert phi'

def sim_pert_theta(P):
    time, Xe, U = default_sim_env(P)
    X0 = np.array(Xe)
    phi0, theta0, psi0 = 0., np.deg2rad(1.), 0.
    X0[fdm.sv_slice_quat] = pal.quat_of_euler([phi0, theta0, psi0])
    # pdb.set_trace()
    return time, sim_open_loop(P, time, X0, U), U, 'pert theta'

def sim_step_up(P, tf=10, dt=0.01):
    time, Xe, U = default_sim_env(P)
    du = np.array([step(t, 0.00005) for t in time])
    U += np.vstack((du, du, -du, -du)).T
    return time, sim_open_loop(P, time, Xe, U), U, 'step up'

def sim_step_uq(P): 
    time, Xe, U = default_sim_env(P)
    du = np.array([step(t, 0.001) for t in time])
    U += np.vstack((du, -du, -du, du)).T
    return time, sim_open_loop(P, time, Xe, U), U, 'step uq'

def sim_step_ur(P):
    time, Xe, U = default_sim_env(P)
    du = np.array([step(t, 0.001) for t in time])
    U += np.vstack((du, -du, du, -du)).T
    return time, sim_open_loop(P, time, Xe, U), U, 'step ur'
        
def sim_step_uz(P):
    time, Xe, U = default_sim_env(P)
    du = np.array([step(t, 0.001) for t in time])
    U += np.vstack((du, du, du, du)).T
    return time, sim_open_loop(P, time, Xe, U), U, 'step uz'

def main():
    np.set_printoptions(linewidth=500)
    P = fdm.Param()
    #test_dyn(P)
    for _s in [sim_pert_phi, sim_pert_theta, sim_step_up, sim_step_uq, sim_step_ur, sim_step_uz]:
    #for _s in [sim_step_up,sim_pert_phi]:
        time, X, U, desc = _s(P)
        fdm.plot(time, X, U, window_title=desc)
    plt.show()
    
    
if __name__ == "__main__":
    main()
