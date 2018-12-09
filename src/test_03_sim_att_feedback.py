#! /usr/bin/env python
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.algebra as pal
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl
import real_time_sim as rts

def main(tf=10., dt=0.005):
    np.set_printoptions(linewidth=500)
    sim = rts.Sim()
    time = np.arange(0, tf, dt)
    X, U = np.zeros((len(time), fdm.sv_size)), np.zeros((len(time), fdm.iv_size))
    Yc = np.zeros((len(time), ctl.iv_size))
    #_in = ctl.CstInput(0, [np.deg2rad(0.1), 0, 0])
    #_in = ctl.SinZInput()
    _in = ctl.RandomInput()
    #_in = ctl.StepEulerInput(pal.e_phi, _a=np.deg2rad(1.), p=4, dt=1)
    #_in = ctl.StepEulerInput(pal.e_theta)
    #_in = ctl.StepEulerInput(pal.e_psi)
    
    X[0] = sim.reset(time[0])
    for i in range(1, len(time)):
        Yc[i-1] = _in.get(time[i-1])
        sim.Yc = Yc[i-1]
        U[i-1], X[i] = sim.run(time[i])
    Yc[-1] = _in.get(time[-1])
    U[-1], unused = sim.run(time[-1])
    sim.fdm.plot(time, X, U)
    sim.ctl.plot(time, Yc, U)
    plt.show()
    
if __name__ == "__main__":
    main()
