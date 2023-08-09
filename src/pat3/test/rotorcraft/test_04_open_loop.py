#! /usr/bin/env python3
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.algebra as pal, pat3.utils as pmu, pat3.frames as pfr
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl
import pat3.vehicles.rotorcraft.multirotor_trajectory as trj


def test1(dt=0.005):
    #_fdm = fdm.SolidFDM()
    _fdm = fdm.UFOFDM()
    #_fdm = fdm.MR_FDM()
    _trj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[1, 0, 0, 0], duration=2.)
    #_trj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[0, 1, 0, 0], duration=2.)
    #_trj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[0, 0, 1, 0], duration=2.)
    #_trj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[0, 0, 0, np.pi], duration=2.)
    #_trj = trj.Circle(c=[0, 0, 0], r=1., v=2., alpha0=0, dalpha=2*np.pi, zt=None, psit=None)
    #_trj = trj.CircleWithIntro(c=[0, 0, 0], r=1., v=3., dt_intro=1.8, dt_stay=0.)
    time = np.arange(0, _trj.duration, dt)
    Yc = np.zeros((len(time), trj._ylen, trj._nder))
    Xr, Ur = np.zeros((len(time), fdm.sv_size)), np.zeros((len(time), fdm.iv_size))
    X = np.zeros((len(time), fdm.sv_size))
    for i in range(0, len(time)):
        Yc[i] = _trj.get(time[i])
        # get reference state ( as euclidian quat) and input ( as Ut, Up, Uq, Ur )
        Xr[i], Ur[i], Xd = ctl.DiffFlatness.state_and_cmd_of_flat_output(None, Yc[i], _fdm.P)
        if i==0:
            X[0] = _fdm.reset(Xr[0], time[0], None) # start at reference state
        else:
            #Usolid = np.array([0, 0, -Ur[i, 0], Ur[i, 1], Ur[i, 2], Ur[i, 3]])
            #X[i] = _fdm.disc_dyn(X[i-1], time[i-1], Usolid, time[i]-time[i-1])

            Uufo = np.array([Ur[i, 0], Ur[i, 1], Ur[i, 2], Ur[i, 3]])
            X[i] = _fdm.disc_dyn(X[i-1], time[i-1], Uufo, time[i]-time[i-1])

            #Xe, Ue = _fdm.trim()
            #X[i] = _fdm.run(None, time[i], -Ue, atm=None)
            #X[i] = _fdm.disc_dyn(_fdm.X, time[i-1], Ue, time[i]-time[i-1])

            
            
        
    #trj.plot(time, Yc)
    fig, axes = fdm.plot(time, Xr, window_title="State Trajectory")#, Ur)
    #pdb.set_trace()
    fdm.plot(time, X, figure=fig, axes=axes)
    
    plt.show()

    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    #main()
    test1()

    
