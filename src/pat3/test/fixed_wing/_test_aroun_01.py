#! /usr/bin/env python3
import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil
import pat3.vehicles.fixed_wing.guidance as p3_guid
import pat3.vehicles.fixed_wing.guidance_ardusoaring as p3_guidardu
import pat3.vehicles.fixed_wing.guidance_soaring as p3_guidsoar


import pat3.atmosphere as p3_atm
import pat3.sensors as p3_sens
import pat3.frames as p3_fr
import pat3.trajectory_3D as p3_traj3d
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.algebra as p3_alg

import pat3.test.fixed_wing.test_03_guidance as t3g

#
# Thermal centering
#
def test_thermal_centering(dm, trim_args, force_recompute, dt=0.01, tf=170.2):
    ref_traj = None
    save_filename = '/tmp/pat_glider_climb.npz'
    atm = p3_atm.AtmosphereWharington(center=[-40., 0, 0], radius=40, strength=-1.5)
    cc = (20., 15., 0.); p0 = 5, 0, -200
    trim_args['va']=9.
    ctl = p3_guidsoar.GuidanceSoaring(dm, ref_traj, trim_args, dt)
    ctl.k_centering = 0.005
    ctl.Xe[0], ctl.Xe[1] , ctl.Xe[2] = p0 # aircraft start point
    ctl.set_circle_center(cc)             # initial circle center
    ctl_logger = p3_guidsoar.Logger()
    sensors = p3_sens.Sensors()
    time, X, U = t3g._run_or_load_sim(dm, ctl, sensors, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename)
    ctl_logger.plot_chronograms(time, X, ctl, atm)
    ctl_logger.plot3D(time, X, ctl, atm)
    ctl_logger.plot_slice_nu(time, X, U, ctl, atm, n0=-80, n1=60, dn=10., h0=180., h1=320., dh=10.)
    return ref_traj, (time, X, U)  



def main(param_filename, trim_args = {'h':0, 'va':11, 'gamma':0}, force_recompute=False):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    ref_traj, (time, X, U) = test_thermal_centering(dm, trim_args, force_recompute)
    plt.show()
    

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
    main(param_filename, force_recompute='-force' in sys.argv)
