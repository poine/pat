#! /usr/bin/env python
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
import pat3.frames as p3_fr
import pat3.trajectory_3D as p3_traj3d
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.algebra as p3_alg

import control.matlab

def run_simulation(dm, ctl, tf=30.5, dt=0.01, trim_args={'h':0, 'va':12, 'gamma':0}, plot=False, atm=None, cbk=None):
    time = np.arange(0, tf, dt)
    X = np.zeros((len(time), dm.sv_size))
    U = np.zeros((len(time),  dm.input_nb()))
    carrots, att_sp = np.zeros((len(time),  3)), np.zeros((len(time), 2))
    X[0] = dm.reset(ctl.Xe, t0=time[0], X_act0=None)#ctl.Ue)
    ctl.enter(X[0], time[0]) # FIXME
    for i in range(1, len(time)):
        Xee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(X[i-1], atm, time[i])
        U[i-1] = ctl.get(time[i-1], X[i-1], Xee)
        carrots[i-1] = ctl.carrot; att_sp[i-1] = ctl.phi_sp, ctl.theta_sp 
        X[i] = dm.run(time[i] - time[i-1], time[i], U[i-1], atm)
        if cbk is not None: cbk()
    Xee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(X[-1], atm, time[-1])
    U[-1] = ctl.get(time[-1], X[-1], Xee); carrots[-1] = ctl.carrot; att_sp[-1] = ctl.phi_sp, ctl.theta_sp
    if cbk is not None: cbk()
    if plot:
        dm.plot_trajectory(time, X, U)
        plt.subplot(5,3,1); plt.plot(time, carrots[:,0])
        plt.subplot(5,3,2); plt.plot(time, carrots[:,1])
        plt.subplot(5,3,3); plt.plot(time, carrots[:,2])
        plt.subplot(5,3,7); plt.plot(time, np.rad2deg(att_sp[:,0]), label='setpoint')
        plt.subplot(5,3,8); plt.plot(time, np.rad2deg(att_sp[:,1]), label='setpoint')
        #plt.show()  
    return time, X, U

def test_line(dm, trim_args, dt=0.01):
    ref_traj = p3_traj3d.LineRefTraj(p1=[0, 10, 0], p2=[400, 10, 0])
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    ctl.v_mode = ctl.v_mode_throttle; ctl.throttle_sp = ctl.Ue[0]
    return ref_traj, run_simulation(dm, ctl, tf=30.5, dt=dt, trim_args=trim_args, plot=True)

def test_circle(dm, trim_args, dt=0.01):
    ref_traj = p3_traj3d.CircleRefTraj(c=[0, 0, 0], r=20)
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    ctl.v_mode = ctl.v_mode_throttle; ctl.throttle_sp = ctl.Ue[0]
    return ref_traj, run_simulation(dm, ctl, tf=30.5, trim_args=trim_args, plot=True)


def test_vctl(dm, trim_args, dt=0.01, compute=False, tf=120.5):
    save_filename = '/tmp/pat_glider_vctl.npz'
    atm = p3_atm.AtmosphereCalm()
    trim_args={'h':26, 'va':17, 'gamma':0}
    #ref_traj = p3_traj3d.CircleRefTraj(c=[70, 0, -25], r=80)
    #ref_traj = p3_traj3d.ZStepCircleRefTraj(c=[70, 0, -30], r=250)
    #ref_traj = p3_traj3d.ZdStepCircleRefTraj(c=[70, 0, -30], r=300)
    ref_traj = p3_traj3d.BankedCircleRefTraj(c=[70, 0, -30], r=150)
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    ctl.v_mode = ctl.v_mode_alt#vz#alt
    ctl_logger = p3_guid.GuidancePurePursuitLogger()
    if compute or not os.path.exists(save_filename):
        time, X, U =  run_simulation(dm, ctl, tf=tf, trim_args=trim_args, plot=False, atm=atm, cbk=lambda:ctl_logger.record(ctl))
        ctl_logger.save(time, X, U, save_filename)
    else:
        time, X, U =  ctl_logger.load(save_filename)

    if 0:
        n0, n1, dn, h0, h1, dh = -20, 100, 5., 0, 60, 2.
        p3_pu.plot_slice_wind_nu(atm, n0=n0, n1=n1, dn=dn, e0=0., h0=h0, h1=h1, dh=dh, zdir=-1.,
                                 show_quiver=True, show_color_bar=True, title="Ridge",
                                 figure=None, ax=None)
        plt.plot(X[:,0],-X[:,2])
    if 1:
        #dm.plot_trajectory_as_ee(time, X, U)
        dm.plot_trajectory(time, X, U)
        Xee = np.array([dm.to_six_dof_euclidian_euler(_X, atm, _t) for _X, _t in zip(X, time)])
        ax=plt.subplot(5,3,3); plt.cla()
        plt.plot(time, -Xee[:,2], label='plant')#dm.sv_zd])
        plt.plot(time, -np.asarray(ctl_logger.carrot)[:,2], label='setpoint')
        p3_pu.decorate(ax, title="$h$", ylab="m", min_yspan=1., legend=True)
        ax=plt.subplot(5,3,6);plt.cla()
        plt.plot(time, -Xee[:,5], label="plant")
        plt.plot(time, -np.asarray(ctl_logger.zd_sp), label="setpoint")
        plt.plot(time, -np.asarray(ctl_logger.zd_sp_traj), label="setpoint_ff")
        p3_pu.decorate(ax, title="$\dot{h}$", ylab="m/s", min_yspan=1., legend=True)
    
    return ref_traj, (time, X, U)  
        
def test_ridge(dm, trim_args, dt=0.01, compute=False):
    save_filename = '/tmp/pat_glider_ridge.npz'
    atm = p3_atm.AtmosphereRidge()
    #atm = p3_atm.AtmosphereCalm()
    trim_args={'h':30, 'va':17, 'gamma':0}
    #ref_traj = p3_traj3d.CircleRefTraj(c=[0, 0, -20], r=20)
    ref_traj = p3_traj3d.BankedCircleRefTraj(c=[120, 0, -40], r=50)
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    #ctl.v_mode = ctl.v_mode_throttle; ctl.throttle_sp = 0.4#ctl.Ue[0]#0.1#
    ctl.v_mode = ctl.v_mode_alt
    ctl_logger = p3_guid.GuidancePurePursuitLogger()
    if compute or not os.path.exists(save_filename):
        time, X, U =  run_simulation(dm, ctl, tf=120.5, trim_args=trim_args, plot=True, atm=atm)
        ctl_logger.save(time, X, U, save_filename)
    else:
        time, X, U =  ctl_logger.load(save_filename)
    #ctl_logger.plot3D(time, X, ctl, atm)
    n0, n1, dn, h0, h1, dh = -20, 100, 5., 0, 60, 2.
    p3_pu.plot_slice_wind_nu(atm, n0=n0, n1=n1, dn=dn, e0=0., h0=h0, h1=h1, dh=dh, zdir=-1.,
                             show_quiver=True, show_color_bar=True, title="Ridge",
                             figure=None, ax=None)
    plt.plot(X[:,0],-X[:,2])
    return ref_traj, (time, X, U)  

def test_shearX(dm, trim_args, dt=0.01, compute=False):
    save_filename = '/tmp/pat_glider_shearX.npz'
    atm = p3_atm.AtmosphereShearX(wind1=15.0, wind2=-2.0, xlayer=60.0, zlayer=40.0)

    trim_args={'h':30, 'va':17, 'gamma':0}

    ref_traj = p3_traj3d.BankedCircleRefTraj(c=[0, 0, -40], r=40)
    ctl = p3_guid.GuidancePurePursuit(dm, ref_traj, trim_args, dt)
    #ctl.v_mode = ctl.v_mode_throttle; ctl.throttle_sp = 0.4#ctl.Ue[0]#0.1#
    ctl.v_mode = ctl.v_mode_alt
    ctl_logger = p3_guid.GuidancePurePursuitLogger()
    if compute or not os.path.exists(save_filename):
        time, X, U =  run_simulation(dm, ctl, tf=35.5, trim_args=trim_args, plot=True, atm=atm)
        ctl_logger.save(time, X, U, save_filename)
    else:
        time, X, U =  ctl_logger.load(save_filename)

    # n0, n1, dn, h0, h1, dh = -20, 100, 5., 0, 60, 2.
    n0, n1, dn, h0, h1, dh = -80, 80., 5., -20, 90, 2.
    p3_pu.plot_slice_wind_nu(atm, n0=n0, n1=n1, dn=dn, e0=0., h0=h0, h1=h1, dh=dh, zdir=-1.,
                             show_quiver=True, show_color_bar=True, title="ShearX",
                             figure=None, ax=None, use_wx=True)
    plt.plot(X[:,0],-X[:,2])
    return ref_traj, (time, X, U)  


def test_climb(dm, trim_args, dt=0.01, compute=True):
    ref_traj = None
    save_filename = '/tmp/pat_glider_climb.npz'
    #atm = p3_atm.AtmosphereWharington(center=[-40., 0, 0], radius=40, strength=-1.5)
    #atm = p3_atm.AtmosphereThermal1()
    #atm.set_params(0., 0., 2000., 300.)
    nc_f = '/home/poine/work/glider_experiments/data/extr_IHODC.1.RK4DI.007.nc'
    atm = p3_atm.AtmosphereNC(nc_f)
    trim_args['va']=9.
    ctl = p3_guidsoar.GuidanceSoaring(dm, ref_traj, trim_args, dt)
    ctl.k_centering = 0.005
    print(ctl.Xe)
    ctl.Xe[0], ctl.Xe[1] , ctl.Xe[2] = 100, 100, -200
    ctl.set_circle_center((50., 15., 0.))
    #ctl.set_circle_center((0., 15., 0.))
    ctl_logger = p3_guidsoar.Logger()
    if compute or not os.path.exists(save_filename):
        time, X, U = run_simulation(dm, ctl, tf=150.02, dt=dt, trim_args=trim_args, plot=True, atm=atm, cbk=lambda:ctl_logger.record(ctl))
        ctl_logger.save(time, X, U, save_filename)
    else:
        time, X, U =  ctl_logger.load(save_filename)
    ctl_logger.plot_chronograms(time, X, ctl, atm)
    ctl_logger.plot3D(time, X, ctl, atm)
    return ref_traj, (time, X, U)  

    #return ref_traj, run_simulation(dm, ctl, tf=120.5, dt=dt, trim_args=trim_args, plot=True, atm=atm, cbk=lambda:ctl_logger.log(ctl))

def test_ardu(dm, trim_args, dt=0.01, compute=False):
    ref_traj = None
    save_filename = '/tmp/pat_glider_ardu.npz'
    ctl_logger = p3_guidardu.Logger()

    #atm = p3_atm.AtmosphereThermal1()
    #atm.set_params(0., 0., 2000., 300.)
    atm = p3_atm.AtmosphereWharington(center=[-40., 0, 0], radius=40, strength=-1.5)
    #atm =  p3_atm.AtmosphereCalm()
    #atm =  p3_atm.AtmosphereCstWind([0, 0, -1.])
    trim_args['va']=9.
    ctl = p3_guidardu.GuidanceArduSoaring(dm, ref_traj, trim_args, dt)

    if compute or not os.path.exists(save_filename):
        time, X, U = run_simulation(dm, ctl, tf=120.02, dt=dt, trim_args=trim_args, plot=True, atm=atm, cbk=lambda:ctl_logger.log(ctl))
        ctl_logger.save(time, X, U, save_filename)
    else:
        time, X, U =  ctl_logger.load(save_filename)
    ctl_logger.plot_chronograms(time, X, ctl, atm)
    ctl_logger.plot2D(time, X, ctl, atm)
    ctl_logger.plot3D(time, X, ctl, atm)
    #Xee = np.array([p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(_X) for _X in X])
    #dm.plot_trajectory_as_ee(time, X, U)
    return ref_traj, (time, X, U)  

def main(param_filename, trim_args = {'h':0, 'va':11, 'gamma':0}, force_recompute=False):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    #ref_traj, (time, X, U) = test_line(dm, trim_args)
    #ref_traj, (time, X, U) = test_circle(dm, trim_args)
    # ref_traj, (time, X, U) = test_vctl(dm, trim_args, compute=force_recompute)
    # ref_traj, (time, X, U) = test_ridge(dm, trim_args, compute=force_recompute)
    #ref_traj, (time, X, U) = test_climb(dm, trim_args, compute=False)
    #ref_traj, (time, X, U) = test_ardu(dm, trim_args)
    ref_traj, (time, X, U) = test_shearX(dm, trim_args, compute=force_recompute)
    # if 1:
    #     p3_pu.plot_3D_traj(ref_traj, X)

    if 0:
        savefile_name = '/tmp/pat_glider_circle.npz'
        np.savez(savefile_name, time=time, X=X)
        print('saved {}'.format(savefile_name))
                
    plt.show()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
    main(param_filename, force_recompute='-force' in sys.argv)
