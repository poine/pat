#! /usr/bin/env python3
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

import pat3.vehicles.fixed_wing.guidance as p3_guid

import control.matlab

import test_03_guidance


def main(param_filename, force_recompute=False, save_filename = '/tmp/pat_glider_wind.npz'):
    
    dm = p1_fw_dyn.DynamicModel(param_filename)
    trim_args = {'h':0, 'va':11, 'gamma':0}
    ref_traj = p3_traj3d.LineRefTraj(p1=[0, 0, 0], p2=[250, 0, 0])
    tf=5.5
    #atm = None
    #atm = p3_atm.AtmosphereCalm()
    #atm = p3_atm.AtmosphereCstWind([-1, 0, 0])
    #atm = p3_atm.AtmosphereVgradient([0, 0, 3])
    #atm = p3_atm.AtmosphereSinetWind([1, 0, 0])
    atm = p3_atm.AtmosphereSinetWind([-5, 0, -1])
    #atm = p3_atm.AtmosphereSinedWind([-5, 0, -1])
    #atm = p3_atm.AtmosphereThermal1()
    dt = 0.01
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    ctl.set_vmode(ctl.v_mode_throttle, ctl.Ue[0])
    #ctl.att_ctl.reset(0, ctl.Xe) # Xe must be SixDOFAeroEuler
    ctl_logger = p3_guid.GuidancePurePursuitLogger()
    #p3_pu.plot_3D_wind(atm)
    #p3_pu.plot_slice_wind(atm)
    #plt.show()
    #time, X, U = test_03_guidance.run_simulation(dm, ctl, None, tf=tf, dt=dt, trim_args=trim_args, atm=atm,
    #                                             cbk=lambda:ctl_logger.record(ctl))
    time, X, U = test_03_guidance._run_or_load_sim(dm, ctl, None, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename)
    ctl_logger.plot_chronograms(time, X, U, ctl, atm)
    #ctl_logger.plot_debug_chronograms(time, X, U, ctl, atm)
    plt.show()

    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = p3_u.pat_ressource('data/vehicles/cularis.xml')
    main(param_filename, force_recompute='-force' in sys.argv)
