#! /usr/bin/env python
import math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

''' Testing fixed wing in open loop '''
import pat3.dynamics as pat_dyn
import pat3.vehicles.fixed_wing.simple_6dof_fdm as fw_dyn

def get_sim_env(tf=0.1, dt=0.01):
    P = fw_dyn.Param(filename='../../../../data/vehicles/cularis.xml')
    #_fdm = fw_dyn.FDM(P)
    _fdm = fw_dyn.FDM3(P)
    Xe, Ue = _fdm.trim({'h':0, 'va':15}, debug=True)
    print(Xe, Ue)
    time = np.arange(0, tf, dt)
    return _fdm, Xe, Ue, time

def test_dyn(_fdm, X0, U, time):
    X = np.zeros((len(time), _fdm.sv_size))
    X[0] = X0
    for i in range(1, len(time)):
        X[i] = _fdm.disc_dyn(X[i-1], time[i-1], U[i-1], time[i]-time[i-1])
    _fdm.plot(time, X)
    plt.show()

def test_equilibrium(tf=0.4, dt=0.01):
    _fdm, Xe, Ue, time = get_sim_env(tf, dt)
    return
    U = np.zeros((len(time), _fdm.iv_size))
    U[:] = Ue
    test_dyn(_fdm, Xe, U, time)

def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    test_equilibrium()
    
    
if __name__ == "__main__":
    main()
