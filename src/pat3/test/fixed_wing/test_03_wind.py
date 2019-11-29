#! /usr/bin/env python
import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.atmosphere as p3_atm
import pat3.trajectory_3D as p3_traj3d
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.algebra as p3_alg

import control.matlab

import test_03_guidance


def main(param_filename):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    trim_args = {'h':0, 'va':11, 'gamma':0}
    ref_traj = p3_traj3d.LineRefTraj(p1=[0, 0, 0], p2=[200, 0, 0])
    #atm = p3_atm.AtmosphereCalm()
    #atm = p3_atm.AtmosphereCstWind([0, 0, 3])
    #atm = p3_atm.AtmosphereVgradient()
    #atm = p3_atm.AtmosphereThermal()
    atm = p3_atm.AtmosphereThermal1()
    
    #p3_pu.plot_3D_wind(atm)
    p3_pu.plot_slice_wind(atm)
    plt.show()
    time, X, U = test_03_guidance.run_simulation(dm, ref_traj, tf=10.5, dt=0.01, trim_args=trim_args,
                                                 plot=True, atm=atm)
    

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
    main(param_filename)
