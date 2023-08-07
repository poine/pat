#! /usr/bin/env python3
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.algebra as pal

''' Testing quadrotor in open loop '''

def test_dyn():
    _fdm = fdm.MR_FDM()
    Xe, Ue = _fdm.trim()
    
    dU = 0.001*np.array([1, -1, -1, 1])
    Xk, Uk = Xe, Ue
    phi0, theta0, psi0 = np.deg2rad(1.), 0., 0.
    Xk[fdm.sv_slice_quat] = pal.quat_of_euler([phi0, theta0, psi0])
    Xkp1 = _fdm.disc_dyn(Xk, 0, Uk, 0.01)
    pdb.set_trace()

def step(t, a=-1., p=10., dt=0.): return a if math.fmod(t+dt, p) > p/2 else -a

def default_sim_env(tf=10, dt=0.01):
    _fdm = fdm.MR_FDM()
    time = np.arange(0, tf, dt) 
    Xe, Ue = _fdm.trim()
    return _fdm, time, Xe, Ue*np.ones((len(time), fdm.iv_size))
    
def sim_open_loop(_fdm, time, X0, U):
    X = np.zeros((len(time), fdm.sv_size))
    X[0] = X0
    for i in range(1, len(time)):
        X[i] = _fdm.disc_dyn(X[i-1], time[i-1], U[i-1], time[i]-time[i-1]) 
    return X

def sim_pert_phi():
    _fdm, time, Xe, U = default_sim_env()
    X0 = np.array(Xe)
    phi0, theta0, psi0 = np.deg2rad(1.), 0., 0.
    X0[fdm.sv_slice_quat] = pal.quat_of_euler([phi0, theta0, psi0])
    # pdb.set_trace()
    return time, sim_open_loop(_fdm, time, X0, U), U, 'pert phi'

def sim_pert_theta():
    _fdm, time, Xe, U = default_sim_env()
    X0 = np.array(Xe)
    phi0, theta0, psi0 = 0., np.deg2rad(1.), 0.
    X0[fdm.sv_slice_quat] = pal.quat_of_euler([phi0, theta0, psi0])
    # pdb.set_trace()
    return time, sim_open_loop(_fdm, time, X0, U), U, 'pert theta'

def sim_step_up(tf=10, dt=0.01):
    _fdm, time, Xe, U = default_sim_env()
    du = np.array([step(t, 0.00005) for t in time])
    U += np.vstack((du, du, -du, -du)).T
    return time, sim_open_loop(_fdm, time, Xe, U), U, 'step up'

def sim_step_uq(): 
    _fdm, time, Xe, U = default_sim_env()
    du = np.array([step(t, 0.00005) for t in time])
    U += np.vstack((du, -du, -du, du)).T
    return time, sim_open_loop(_fdm, time, Xe, U), U, 'step uq'

def sim_step_ur():
    _fdm, time, Xe, U = default_sim_env()
    du = np.array([step(t, 0.001) for t in time])
    U += np.vstack((du, -du, du, -du)).T
    return time, sim_open_loop(_fdm, time, Xe, U), U, 'step ur'
        
def sim_step_uz():
    _fdm, time, Xe, U = default_sim_env()
    du = np.array([step(t, 0.001) for t in time])
    U += np.vstack((du, du, du, du)).T
    return time, sim_open_loop(_fdm, time, Xe, U), U, 'step uz'

def main():
    np.set_printoptions(linewidth=500)
    if 0:
        test_dyn()
    if 1:
        for _s in [sim_pert_phi, sim_pert_theta, sim_step_up, sim_step_uq, sim_step_ur, sim_step_uz]:
        #for _s in [sim_step_up,sim_pert_phi]:
            time, X, U, desc = _s()
            fdm.plot(time, X, U, window_title=desc)
        plt.show()
    
    
if __name__ == "__main__":
    main()
