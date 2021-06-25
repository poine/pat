#! /usr/bin/env python3
'''
Unit test for fixed wing guidance
'''
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

import control.matlab

LOG = logging.getLogger(__name__)

def run_simulation(dm, ctl, sensors=None, tf=30.5, dt=0.01, trim_args={'h':0, 'va':12, 'gamma':0}, atm=None, cbk=None, desc=None):
    report = f'\n  Running simulation for {tf}s'
    if desc is not None: report += '\n  desc: '+desc
    LOG.info(report)
    time = np.arange(0, tf, dt)
    X, U = np.zeros((len(time), dm.sv_size)), np.zeros((len(time),  dm.input_nb()))
    X[0] = dm.reset(ctl.Xe, t0=time[0], X_act0=None)#ctl.Ue)
    ctl.enter(X[0], time[0]) # FIXME: fix what?
    for i in range(1, len(time)):
        Xee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(X[i-1], atm, time[i])
        if sensors is None: Xaes, Xees = X[i-1], Xee
        else: Xaes, Xees = sensors.get_measurements(X[i-1], Xee)
        U[i-1] = ctl.get(time[i-1], Xaes, Xees)
        if cbk is not None: cbk()  # used for recording ctl variables
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1], atm)
    Xee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(X[-1], atm, time[-1])
    U[-1] = ctl.get(time[-1], X[-1], Xee)
    if cbk is not None: cbk()
    return time, X, U

def test_line(dm, trim_args, dt=0.01):
    ref_traj = p3_traj3d.LineRefTraj(p1=[0, 10, 0], p2=[400, 10, 0])
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    ctl.v_mode = ctl.v_mode_throttle; ctl.throttle_sp = ctl.Ue[0]
    return ref_traj, run_simulation(dm, ctl, tf=30.5, dt=dt, trim_args=trim_args)

def test_circle(dm, trim_args, dt=0.01):
    ref_traj = p3_traj3d.CircleRefTraj(c=[0, 0, 0], r=20)
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    ctl.v_mode = ctl.v_mode_throttle; ctl.throttle_sp = ctl.Ue[0]
    return ref_traj, run_simulation(dm, ctl, tf=30.5, trim_args=trim_args)

# load and save simulations
def _run_or_load_sim(dm, ctl, sensors, tf, dt, trim_args, atm, ctl_logger, force_compute, save_filename, desc=None):
    if force_compute or not os.path.exists(save_filename):
        time, X, U =  run_simulation(dm, ctl, sensors, tf=tf, dt=dt, trim_args=trim_args, atm=atm, cbk=lambda:ctl_logger.record(ctl), desc=desc)
        ctl_logger.save(time, X, U, save_filename)
    else:
        time, X, U =  ctl_logger.load(save_filename)
    return time, X, U

#
# test vertical control on a few trajectories
#
def test_vctl(dm, trim_args, force_recompute=False, dt=0.01, tf=40.5):
    save_filename = '/tmp/pat_glider_vctl_w.npz'
    #atm = p3_atm.AtmosphereCalm()
    atm = p3_atm.AtmosphereCstWind([1, 0, 0])
    trim_args={'h':26, 'va':17, 'gamma':0}
    ref_traj = p3_traj3d.CircleRefTraj(c=[70, 0, -25], r=80)
    #ref_traj = p3_traj3d.ZStepCircleRefTraj(c=[70, 0, -30], r=250)
    #ref_traj = p3_traj3d.ZdStepCircleRefTraj(c=[70, 0, -30], r=300)
    #ref_traj = p3_traj3d.BankedCircleRefTraj(c=[70, 0, -30], r=150)
    #ref_traj = p3_traj3d.SquareRefTraj()
    #ref_traj = p3_traj3d.OctogonRefTraj() # kindof broken
    #ref_traj = p3_traj3d.OvalRefTraj() #
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    ctl.v_mode = ctl.v_mode_alt#vz#alt
    ctl_logger = p3_guid.GuidancePurePursuitLogger()

    time, X, U = _run_or_load_sim(dm, ctl, None, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename)

    ctl_logger.plot3D(time, X, ctl, atm, ref_traj)
    ctl_logger.plot_slice_nu(time, X, U, ctl, atm)
    ctl_logger.plot_chronograms(time, X, U, ctl, atm)
        
    return ref_traj, (time, X, U)  

#
# Slope soaring, gliding with throttle off
#
def test_slope_soaring(dm, trim_args, force_recompute=False, dt=0.01, tf=120.5):
    save_filename = '/tmp/pat_glider_slope_soaring.npz'
    atm = p3_atm.AtmosphereRidge()
    trim_args={'h':30, 'va':12, 'gamma':0}
    #ref_traj = p3_traj3d.CircleRefTraj(c=[-10, 0, -50], r=25)
    ref_traj = p3_traj3d.OvalRefTraj( _r=25., _l=200, _c=np.array([-10, 0, -50])) #
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    ctl.v_mode = ctl.v_mode_throttle; ctl.throttle_sp = 0. # we glide
    ctl_logger = p3_guid.GuidancePurePursuitLogger()

    time, X, U = _run_or_load_sim(dm, ctl, None, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename)

    ctl_logger.plot3D(time, X, ctl, atm, ref_traj)
    ctl_logger.plot_slice_nu(time, X, U, ctl, atm, n0=-50, n1=40)
    ctl_logger.plot_chronograms(time, X, U, ctl, atm)

    return ref_traj, (time, X, U)

#
# Dynamic soaring
#

def _ds_exp(i):
    _exp = [(p3_atm.AtmosphereCalm(), '/tmp/pat_glider_ds_nw.npz'),
            (p3_atm.AtmosphereCstWind([1, 0, 0]), '/tmp/pat_glider_ds_wc100.npz'),
            (p3_atm.AtmosphereRidge(), '/tmp/pat_glider_ds_wr.npz'),
            (p3_atm.AtmosphereShearX(wind1=5.0, wind2=0.0, xlayer=60.0, zlayer=40.0), '/tmp/pat_glider_ds_ws50.npz'),
            (p3_atm.AtmosphereVgradient(w0=-2, w1=5, h0=20, h1=60), '/tmp/pat_glider_ds_wvg.npz'),
            (p3_atm.AtmosphereVgradient(w0=-7.5, w1=7.5, h0=20, h1=60), '/tmp/pat_glider_ds_wvg2.npz')]
    atm, save_filename = _exp[i]
    return atm, save_filename
        
def test_dynamic_soaring(dm, trim_args, force_recompute=False, dt=0.005, tf=150.):
    #atm = p3_atm.AtmosphereShearX(wind1=15.0, wind2=-2.0, xlayer=60.0, zlayer=40.0)
    #atm = p3_atm.AtmosphereShearX(wind1=7.0, wind2=-1.0, xlayer=60.0, zlayer=40.0)
    #atm = p3_atm.AtmosphereShearX(wind1=5.0, wind2=0.0, xlayer=60.0, zlayer=40.0)
    atm, save_filename = _ds_exp(4)
    trim_args={'h':30, 'va':17, 'gamma':0}
    #ref_traj = p3_traj3d.CircleRefTraj(c=[0, 0, -20], r=20)
    ref_traj = p3_traj3d.BankedCircleRefTraj(c=[100, 0, -40], r=60, slope=np.deg2rad(10))
    #ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt, lookahead_dist=15., max_phi=np.deg2rad(45))
    ctl = p3_guid.GuidanceDS(dm, ref_traj, trim_args, dt, lookahead_dist=15., max_phi=np.deg2rad(60))
    ctl.v_mode = ctl.v_mode_alt
    ctl_logger = p3_guid.GuidancePurePursuitLogger()
    time, X, U = _run_or_load_sim(dm, ctl, None, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename)

    ctl_logger.plot3D(time, X, ctl, atm)
    ctl_logger.plot_slice_nu(time, X, U, ctl, atm, ref_traj, n0=0, n1=210, dn=10, h0=0, h1=80, dh=10)
    ctl_logger.plot_chronograms(time, X, U, ctl, atm)
    return ref_traj, (time, X, U)  

#
# Thermal centering
#
def test_thermal_centering(dm, trim_args, force_recompute, dt=0.01, tf=170.2):
    ref_traj = None
    exp=1
    save_filename = f'/tmp/pat_glider_climb_{exp}.npz'
    if exp==0:
        atm = p3_atm.AtmosphereWharington(center=[-40., 0, 0], radius=40, strength=-1.5)
        cc = (10., 15., 0.); p0 = 5, 0, -200
        n_spec, e0, h_spec = (-100, 50, 10), 0., (180., 350., 20.)
    elif exp==1:
        nc_f = '/home/poine/work/glider_experiments/data/extr_IHODC.1.RK4DI.007.nc'
        atm = p3_atm.AtmosphereNC(nc_f)
        cc = (80., 100., -100.); p0 = 70, 100, -100
        n_spec, e0, h_spec = (0, 400, 10), 100., (100., 350., 20.)
    trim_args['va']=9.
    ctl = p3_guidsoar.GuidanceSoaring(dm, ref_traj, trim_args, dt)
    ctl.k_centering = 0.015
    ctl.Xe[0], ctl.Xe[1] , ctl.Xe[2] = p0 # aircraft start point
    ctl.set_circle_center(cc)             # initial circle center
    ctl_logger = p3_guidsoar.Logger()
    sensors = p3_sens.Sensors(std_va=0.1, std_vz=0.05)
    time, X, U = _run_or_load_sim(dm, ctl, sensors, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename)
    ctl_logger.plot_chronograms(time, X, ctl, atm)
    ctl_logger.plot3D(time, X, ctl, atm)
    ctl_logger.plot_slice_nu(time, X, U, ctl, atm, n_spec, e0, h_spec)
    return ref_traj, (time, X, U)  


#
# Ardusoaring Thermal centering
#
def test_ardusoaring(dm, trim_args, force_recompute, dt=0.01, tf=170.2):
    ref_traj = None
    exp=1
    save_filename = f'/tmp/pat_glider_ardu_{exp}.npz'
    if exp==0:
        atm = p3_atm.AtmosphereWharington(center=[-40., 0, 0], radius=40, strength=-1.5)
        p0 = 5, 0, -100
        n_spec, e0, h_spec = (-100, 50, 10), 0, (180., 350., 20.)
    elif exp==1:
        nc_f = '/home/poine/work/glider_experiments/data/extr_IHODC.1.RK4DI.007.nc'
        atm = p3_atm.AtmosphereNC(nc_f)
        p0 = 70, 100, -100
        n_spec, e0, h_spec = (0, 400, 10), 100., (100., 350., 20.)
    trim_args['va']=9.
    ctl = p3_guidardu.GuidanceArduSoaring(dm, ref_traj, trim_args, dt)
    #ctl.reset(np.pi)#, xc0=None, yc0=None, s0=INITIAL_THERMAL_STRENGTH, r0=INITIAL_THERMAL_RADIUS)
    ctl.Xe[0], ctl.Xe[1] , ctl.Xe[2] = p0 # aircraft start point
    sensors = p3_sens.Sensors(std_va=0.1, std_vz=0.05)
    ctl_logger = p3_guidardu.Logger()
    time, X, U = _run_or_load_sim(dm, ctl, sensors, tf, dt, trim_args, atm, ctl_logger, force_recompute, save_filename)
 
    ctl_logger.plot_chronograms(time, X, ctl, atm)
    ctl_logger.plot2D(time, X, ctl, atm, n_spec, e0, h_spec)
    ctl_logger.plot3D(time, X, ctl, atm)
    return ref_traj, (time, X, U)  

def main(param_filename, trim_args = {'h':0, 'va':11, 'gamma':0}, force_recompute=False):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    #ref_traj, (time, X, U) = test_line(dm, trim_args)
    #ref_traj, (time, X, U) = test_circle(dm, trim_args)
    #ref_traj, (time, X, U) = test_vctl(dm, trim_args, force_recompute)
    #ref_traj, (time, X, U) = test_slope_soaring(dm, trim_args, force_recompute)
    #ref_traj, (time, X, U) = test_dynamic_soaring(dm, trim_args, force_recompute)
    #ref_traj, (time, X, U) = test_thermal_centering(dm, trim_args, force_recompute)
    ref_traj, (time, X, U) = test_ardusoaring(dm, trim_args, force_recompute)
    if 0: p3_pu.plot_3D_traj(ref_traj, X)

    plt.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
    main(param_filename, force_recompute='-force' in sys.argv)
