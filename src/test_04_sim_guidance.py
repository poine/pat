#! /usr/bin/env python

import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.algebra as pal, pat3.utils as pmu
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl
import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt

def main(traj, dt=0.005):
    _fdm = fdm.FDM()
    _ctl_input = ctl.TrajRef(traj, _fdm)
    _ctl = ctl.PosController(_fdm, _ctl_input)
    sim = pmu.Sim(_fdm, _ctl)

    time = np.arange(0, traj.duration, dt)
    X, U = np.zeros((len(time), fdm.sv_size)), np.zeros((len(time), fdm.iv_size))
    Xref = np.zeros((len(time), _ctl.ref_size))
    X[0] = sim.reset(time[0], _ctl_input.get(time[0])[2])
    for i in range(1, len(time)):
        U[i-1], X[i] = sim.run(time[i])
        Xref[i-1] = _ctl.Xref
    U[-1], Xref[-1] = U[-2], Xref[-2]
        
    figure = _fdm.plot(time, X, U)
    _ctl.plot(time, U, Xref)
    plt.show()
    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    #_traj =  pmt.Circle(v=1.)
    _traj =  pmt.FigureOfHeight(v=1.)
    main(_traj)
