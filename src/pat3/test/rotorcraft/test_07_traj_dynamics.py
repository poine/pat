#! /usr/bin/env python3
'''
  Testing space indexed trajectory
'''
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.vehicles.rotorcraft.multirotor_trajectory as trj
import pat3.vehicles.rotorcraft.multirotor_trajectory_dev as trj_dev


import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl

import pdb

def plot_dynamics(_dyn, time, fig=None, axs=None, window_title="Trajectory Dynamics"):
    lambdas = np.array([_dyn.get(t) for t in time])
    if fig is None: fig, axs = plt.subplots(trj._nder, 1)
    if window_title is not None: fig.canvas.manager.set_window_title(window_title)
    for i in range(trj._nder):
        axs[i].plot(time, lambdas[:,i])
    return fig, axs

def main(duration=10., dt=1./200):
    _fdm = fdm.MR_FDM()
    straj = trj_dev.SpaceCircle(r=1.5, c=[0,1.5], alpha0=0, dalpha=2*np.pi)
    _f1, _a1, _f2, _a2, _f3, _a3 = [None]*6
    for duration in [4, 6, 8, 10]:
        time = np.arange(0, duration, dt)
        #dtraj = trj.AffineOne(1./duration,0., duration)
        dtraj = trj.PolynomialOne([0,0,0,0,0], [1,0,0,0,0], duration)
        #dtraj = trj.SigmoidOne(duration, 0.5)
        _f1, _a1 = plot_dynamics(dtraj, time, _f1, _a1)
        _trj = trj_dev.SpaceIndexedTraj(straj, dtraj)
        Yc = np.zeros((len(time), trj._ylen, trj._nder))
        for i in range(0, len(time)):
            Yc[i] = _trj.get(time[i])
        _f2, _a2 = trj.plot(time, Yc, figure=_f2, axes=_a2)
        Xr, Ur = np.zeros((len(time), fdm.sv_size)), np.zeros((len(time), fdm.iv_size))
        for i in range(0, len(time)):
            Xr[i], Ur[i], Xd = ctl.DiffFlatness.state_and_cmd_of_flat_output(None, Yc[i], _fdm.P)
        _f3, _a3 = fdm.plot(time, Xr, window_title="State Trajectory", U=None, figure=_f3, axes=_a3)
    plt.show()
  

if __name__ == "__main__":
    main()
