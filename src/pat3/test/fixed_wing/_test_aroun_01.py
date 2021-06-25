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
def get_defaults():
    param_filename =  p3_u.pat_ressource('data/vehicles/cularis.xml')
    dm = p1_fw_dyn.DynamicModel(param_filename)
    atm = p3_atm.AtmosphereWharington(center=[-40., 0, 0], radius=40, strength=-1.5)
    sensors = p3_sens.Sensors(std_va=0.)
    cc = (20., 15., 0.); p0 = 5, 0, -200
    trim_args = {'h':0, 'va':9, 'gamma':0}
    dt=0.01
    ctl = p3_guidsoar.GuidanceSoaring(dm, None, trim_args, dt)
    ctl.Xe[0], ctl.Xe[1] , ctl.Xe[2] = p0 # aircraft start point
    ctl.set_circle_center(cc)             # initial circle center
    return dm, ctl, atm, sensors, dt, p0, cc, trim_args


def test_thermal_centering(defaults, force_recompute, tf=170.2,
                           save_filename='/tmp/pat_glider_climb.npz'):
    dm, ctl, atm, sensors, dt, p0, cc, trim_args = defaults
    ctl.k_centering = 0.005#0.01#0.005
    ctl_logger = p3_guidsoar.Logger()
    time, X, U = t3g._run_or_load_sim(dm, ctl, sensors, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename)
    ctl_logger.plot_chronograms(time, X, ctl, atm)
    f, ax = ctl_logger.plot3D(time, X, ctl, atm)
    #plt.show() # prints the 3D view orientation after you set it manually
    #print(f'ax.azim {ax.azim} ax.elev {ax.elev}')
    ctl_logger.plot_slice_nu(time, X, U, ctl, atm, n0=-80, n1=60, dn=10., h0=180., h1=320., dh=10.)
    return time, X, U


def get_figs_and_axes(n_cases):
    fig1 = plt.figure(tight_layout=True, figsize=[4.*n_cases, 3*2.])
    axes1 = fig1.subplots(3, n_cases, sharex=True).T
    fig2 = plt.figure(tight_layout=True, figsize=[4.*n_cases, 4.])
    axes2 = fig2.subplots(1, n_cases, sharey=True)
    return fig1, axes1, fig2, axes2


def test_noises(defaults, force_recompute, tf=170.2,
                save_fmt='/tmp/pat_glider_climb_noise_{:.1e}.npz'):
    dm, ctl, atm, sensors, dt, p0, cc, trim_args = defaults
    noises, ctl.k_centering = [0., 0.025, 0.05, 0.1], 0.015# 0.006
    fig1, axes1, fig2, axes2 = get_figs_and_axes(len(noises))
    for i, (noise, ax1, ax2) in enumerate(zip(noises, axes1, axes2)):
        sensors.vz_std = noise
        save_filename = save_fmt.format(noise)
        ctl_logger = p3_guidsoar.Logger()
        desc = f'soaring noise {noise:.2f}m/s gain {ctl.k_centering:.1e}'
        ctl.set_circle_center(cc)             # initial circle center... that sucks :(
        time, X, U = t3g._run_or_load_sim(dm, ctl, sensors, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename, desc)
        title = f'noise: {noise:f}'
        ctl_logger.plot_chronograms(time, X, ctl, atm, title, None, fig1, ax1)
        ctl_logger.plot_slice_nu(time, X, U, ctl, atm, n0=-80, n1=60, dn=10., h0=180., h1=320., dh=10., title=title, fig=fig2, ax=ax2)


def test_gains(defaults, force_recompute, tf=170.2,
               save_fmt='/tmp/pat_glider_climb_gain_{:.1e}.npz'):
    dm, ctl, atm, sensors, dt, p0, cc, trim_args = defaults
    gains = [0.0050, 0.015, 0.018, 0.02, 0.03]
    fig1, axes1, fig2, axes2 = get_figs_and_axes(len(gains))
    for i, (gain, ax1, ax2) in enumerate(zip(gains, axes1, axes2)):
        ctl.k_centering = gain
        save_filename, ctl_logger = save_fmt.format(gain), p3_guidsoar.Logger()
        desc = f'soaring noise: {sensors.va_std:.1f} gain {gain:.2e}'
        ctl.set_circle_center(cc)             # initial circle center... that sucks :(
        time, X, U = t3g._run_or_load_sim(dm, ctl, sensors, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename, desc)
        title = f'gain: {gain:f}'
        ctl_logger.plot_chronograms(time, X, ctl, atm, title, None, fig1, ax1)
        ctl_logger.plot_slice_nu(time, X, U, ctl, atm, n0=-80, n1=60, dn=10., h0=180., h1=320., dh=10., title=title, fig=fig2, ax=ax2)



        
def main(force_recompute=False):
    dm, ctl, atm, sensors, dt, p0, cc, trim_args = defaults = get_defaults()
    #time, X, U = test_thermal_centering(defaults, force_recompute)#, tf=50.)
    test_gains(defaults, force_recompute)#, tf=50.)
    #test_noises(defaults, force_recompute)#, tf=50.)
    plt.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    main(force_recompute='-force' in sys.argv)
