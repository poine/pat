#! /usr/bin/env python
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.algebra as pal
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl

def step(t, a=-1., p=10., dt=0.): return a if math.fmod(t+dt, p) > p/2 else -a


def step_z(time): return np.array([[step(t), 1, 0, 0, 0] for t in time])


def step_phi(time):
    phis = [step(t, a=np.deg2rad(1.)) for t in time]
    zc = np.zeros((len(time), 1))
    qc = np.array([pal.quat_of_euler([phi, 0, 0]) for phi in phis])
    Yc = np.hstack((zc, qc))
    return Yc

def step_euler(time, _i, _a=np.deg2rad(1.)):
    eulers = np.zeros((len(time), 3))
    eulers[:,_i] = [step(t, a=np.deg2rad(_a)) for t in time]
    zc, qc = np.zeros((len(time), 1)), [pal.quat_of_euler(_e) for _e in eulers]
    return np.hstack((zc, qc))

def step_euler2(time, _i, _a):
    _in = ctl.StepEulerInput(_i, _a)
    Yc = [_in.get(t) for t in time]
    return np.array(Yc)
    
def plot_sp(time, Yc):
    plt.subplot(5, 3, 3)
    plt.plot(time, Yc[:,0])
    euler_c = np.array([pal.euler_of_quat(_sp[1:]) for _sp in Yc])
    plt.subplot(5, 3, 7)
    plt.plot(time, np.rad2deg(euler_c[:,0]))
    plt.subplot(5, 3, 8)
    plt.plot(time, np.rad2deg(euler_c[:,1]))
    plt.subplot(5, 3, 9)
    plt.plot(time, np.rad2deg(euler_c[:,2]))

def main(tf=10., dt=0.01):
    np.set_printoptions(linewidth=500)
    _fdm = fdm.MR_FDM()
    _ctl = ctl.ZAttController(_fdm)
    time = np.arange(0, tf, dt) 

    #Yc = step_phi(time)
    #Yc = step_euler(time, pal.e_phi, 1.)
    Yc = step_euler(time, pal.e_theta, 1.)
    #Yc = step_euler2(time, pal.e_psi, np.deg2rad(1.))
    #Yc = step_z(time)

    Xe, Ue = _fdm.trim()    
    X0 = np.array(Xe)
    #X0[fdm.sv_zd] += 0.1
    #X0[fdm.sv_slice_quat] = pal.quat_of_euler([np.deg2rad(1.), 0., 0.])
    X, U = np.zeros((len(time), fdm.sv_size)), Ue*np.ones((len(time), fdm.iv_size))
    X[0] = _fdm.reset(X0, time[0])
    for i in range(1, len(time)):
        U[i-1] = _ctl.get(time[i-1], X[i-1], Yc[i-1])
        #X[i] = _fdm.run(time[i], U[i-1])
        X[i] = _fdm.disc_dyn(X[i-1], time[i-1], U[i-1], time[i]-time[i-1])
        
    fdm.plot(time, X, U)
    plot_sp(time, Yc)
    plt.show()
    
if __name__ == "__main__":
    main()
