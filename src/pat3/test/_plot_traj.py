#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#

import numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.utils as p3_u, pat3.atmosphere as p3_atm, pat3.trajectory_3D as p3_traj3d, pat3.frames as p3_fr
import pat3.plot_utils as p3_plu
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.guidance as p3_guid

def main(dt=0.005):
    param_filename = p3_u.pat_ressource('data/vehicles/cularis.xml')
    dm = p1_fw_dyn.DynamicModel(param_filename)
    trim_args={'h':30, 'va':17, 'gamma':0}
    if 0:
        atm = p3_atm.AtmosphereCalm()
        save_filename = '/tmp/pat_glider_ds_nw.npz'
    if 0:
        atm = p3_atm.AtmosphereCstWind([1, 0, 0])
        save_filename = '/tmp/pat_glider_ds_wc100.npz'
    if 0:
        atm = p3_atm.AtmosphereRidge()
        save_filename = '/tmp/pat_glider_ds_wr.npz'
    if 0:
        atm = p3_atm.AtmosphereShearX(wind1=5.0, wind2=0.0, xlayer=60.0, zlayer=40.0)
        save_filename = '/tmp/pat_glider_ds_ws50.npz'
    if 1:
        atm = p3_atm.AtmosphereVgradient()
        save_filename = '/tmp/pat_glider_ds_wvg.npz'

    ref_traj = p3_traj3d.BankedCircleRefTraj(c=[100, 0, -40], r=60, slope=np.deg2rad(10))
    ctl = p3_guid.GuidanceDS(dm, ref_traj, trim_args, dt, lookahead_dist=15., max_phi=np.deg2rad(60))
    ctl_logger = p3_guid.GuidancePurePursuitLogger()
    time, Xae, U =  ctl_logger.load(save_filename)

    # aliasing state variables
    _pos   = Xae[:,p3_fr.SixDOFAeroEuler.sv_slice_pos]
    _alt   = -Xae[:, p3_fr.SixDOFAeroEuler.sv_z]
    _va    = Xae[:, p3_fr.SixDOFAeroEuler.sv_va]
    _va2   = _va**2
    _alpha = Xae[:, p3_fr.SixDOFAeroEuler.sv_alpha]
    _eul  = Xae[:,p3_fr.SixDOFAeroEuler.sv_slice_eul]
    _phi, _theta, _psi = [_eul[:,_i] for _i in range(3)]
    _gamma = _theta - _alpha
    groundspeed_3d = np.linalg.norm(_vi, axis=1) # I hate this variable name
    # aliasing input variables
    _throttle = U[:,dm.iv_dth()]

    # compute state variable in euclidian/euler format
    Xee = np.array([dm.to_six_dof_euclidian_euler(_X, atm, 0.) for _X in Xae])
    _vi = Xee[:, p3_fr.SixDOFEuclidianEuler.sv_slice_vel]
    
    #ctl_logger.plot_chronograms(time, X, U, ctl, atm2)

    # Compute Energy
    e_kin_air = 0.5*dm.P.m*_va2
    alt_0 = 30.
    e_pot = dm.P.m*9.81*(_alt-alt_0)
    e_tot_air = e_kin_air + e_pot

    # Compute wind in body frame
    #df.wx.iloc[i] = -V_a * cos(theta-np.deg2rad(alpha)) + V_g
    #df.wz.iloc[i] = V_a * sin(theta-np.deg2rad(alpha)) + V_z
    #Vg is gps ground speed, in X-Y
    #and Vz is climb speed
    # def estimate_inplane_wind(df_in):
    #     df = df_in.copy() # Be careful, this is important !!!
    #     for i in range(df.index.shape[0]):
    #         V_a   = df.airspeed.iloc[i]
    #         theta = df.theta.iloc[i]
    #         #alpha = df.alpha.iloc[i]
    #         alpha = alpha_func(df.index[i])
    #         V_g = df.vel.iloc[i]
    #         #         V_g = df.vel_3d.iloc[i]
    #         V_z = -df.climb.iloc[i] #
    #         df.wx.iloc[i] = -V_a * cos(theta-np.deg2rad(alpha)) + V_g
    #         df.wz.iloc[i] = V_a * sin(theta-np.deg2rad(alpha)) + V_z # going up is negative for wz
    #         df.alpha[i] = alpha
    #     return df
    Wned = np.array([atm.get_wind_ned(_p, _t) for _p, _t in zip(_pos, time)])
    Wb = np.array([p3_fr.vel_world_to_body_eul(_v, _e) for _v, _e in zip(Wned, _eul)])
    Wx, Wy, Wz = [Wb[:,_i] for _i in range(3)]
    dt = time[1]-time[0]
    dWx, dWz = np.gradient(Wx, dt), np.gradient(Wz, dt) 
    # compute power
    rho = 1.225
    AR = 11.84
    e = 0.85
    # I don't have accel, let's compute it
    Az = np.gradient(_vi[:,2], dt)
    Lift = -Az * dm.P.m
    #Lift = -Az * mass
    CL = Lift/ (0.5*rho*_va2*dm.P.Sref)
    CD = dm.P.CD0 + CL**2/(np.pi*AR*e)
    D = -0.5*rho*_va2*dm.P.Sref*CD
    P_drag= _va * D
    P_dwx = -dWx*(-_va * np.sign(_gamma) * np.cos(_gamma) )
    P_dwz =  dWz*(_va * np.sin(_gamma) )
     
    fig = plt.figure(figsize=(10,14.5))
    axes = fig.subplots(6, 1, sharex=True)
    ax = axes[0]#fig.add_subplot(611)
    ax.set_title('Energy in Air-Path Frame')
    ax.plot(time,e_tot_air, label='$E_{Total}$')
    ax.plot(time,e_kin_air, label='$E_{Kinetic}$')
    ax.plot(time,e_pot, label='$E_{Potential}$')
    ax.grid();ax.legend();ax.set_ylabel('Energy [J]')#plt.ylabel('Energy/(mg) [m]')
    #ax.set_xlim([xlim1,xlim2]);#;plt.xlim([0,1])#;plt.ylim([0,70])
    ax.set_xticklabels([])

    
    ax = axes[1]#fig.add_subplot(612)
    ax.plot(time, _alt-alt_0, label='Height AGL')
    ax.plot(time, _va, label='$V_a$');
    ax.plot(time, groundspeed_3d, label='$V_{i}$')
    ax.grid();
    #ax.set_ylabel('Height AGL [m] $\&$ \\\ Speed [m/s]')
    ax.set_ylabel('Height AGL [m] \\\ Speed [m/s]')
    ax.legend();
    #ax.set_xlim([xlim1,xlim2]); #plt.ylim([0,35]) #plt.xlim([0,1]);


    ax = axes[2]#fig.add_subplot(613)
    ax.plot(time, np.rad2deg(_gamma), label='$\gamma$ [deg]');
    ax.plot(time, Wx, label='$W_x$');
    ax.plot(time,dWx, label='$\dot{W_x}$');
    ax.plot(time,Wz, label='$W_z$');
    ax.plot(time,-_vi[:,2], label='$Vi_z$');

    ax.grid();ax.legend();ax.set_ylabel('  Flight Path ($\gamma$) [deg] \\\  Wind Speed [$m/s$] \\\  Gradient [$m/s^2$]');#plt.ylim([-12,12]);
    #ax.set_xlim([xlim1,xlim2]);ax.set_xticklabels([])


    ax = axes[3]#fig.add_subplot(614)
    #ax.plot(time, dWx, label='dwx'); ax.plot(time, dWz, label='dwz')
    ax.plot(time, P_drag, color='red', label='$P_D$  [W], $\sum{P_D}$ = %0.2f [Ws]' % (np.nansum(P_drag)/100) ) 
    ax.plot(time,P_dwx, color='black', alpha=0.6, label='$P_{\dot{W}_X}$ [W], $\sum{P_{\dot{W}_X}}$ = %0.2f [Ws]' % (np.nansum(P_dwx)/100) );
    ax.fill_between(time,0,P_dwx, where=(P_dwx >= 0), alpha=0.50, color='green', interpolate=True)
    ax.fill_between(time,0,P_dwx, where=(P_dwx < 0), alpha=0.50, color='red', interpolate=True)
    ax.grid();ax.legend()#;ax.set_ylabel('Power ($P_{\dot{W}_X}\, \& \, P_D$)  [W]')#;ax.set_xlim([xlim1,xlim2]);plt.ylim([-100,120])
    ax.set_xticklabels([])

    ax = axes[4]#fig.add_subplot(615)
    ax.plot(time,P_dwz, color='black', alpha=0.3, label='$P_{\dot{W_Z}}$, $\sum{P_{\dot{W_Z}}}$ = %0.2f [Ws]' % (np.nansum(P_dwz)/100) );
    ax.fill_between(time,0,P_dwz, where=(P_dwz >= 0), alpha=0.50, color='green', interpolate=True)#, label='P$_{\dot{w}}$');
    ax.fill_between(time,0,P_dwz, where=(P_dwz < 0), alpha=0.50, color='red', interpolate=True)#, label='P$_{\dot{w}}$');
    ax.grid();plt.legend();plt.ylabel('Power ($P_{\dot{W}_Z}$)  [W]')#;ax.set_xlim([xlim1,xlim2]);ax.set_ylim([-100,120])#;plt.xlim([0,1])#;plt.ylim([-10,10])
    ax.set_xticklabels([])

    ax = axes[5]#fig.add_subplot(616);
    ax.plot(time, _throttle, label='Throttle,  Averaged = %0.2f' % (np.nanmean(_throttle)) );
    #ax.fill_between(time,0,electrical_power, where=(_throttle >= throttle_limit), alpha=0.20, color='grey', interpolate=True)
    #ax.plot(time,electrical_power, label='$P_{Elec.}$ [W] , $\sum{P_{Elec.}}$ = %0.2f [Ws]' % (np.nansum(electrical_power[throttle>=throttle_limit])/100.* propulsion_eff) ); #np.nanmean(electrical_power),
    ax.grid();plt.legend();ax.set_ylabel('Throttle [\%] \& \\\ Elec. Power ($P_{Elec.}$) [W]')#;ax.set_xlim([xlim1,xlim2]);#;plt.xlim([xlim1,xlim2]);plt.ylim([0,100])
    ax.set_xlabel('Time [s]')
    ax.set_xlim([20, 30])
    
    plt.savefig('/tmp/traj_murat_ds.png', dpi=120, bbox_inches='tight')

def plot_atm():
    atm1 =  p3_atm.AtmosphereShearX(wind1=5.0, wind2=0.0, xlayer=60.0, zlayer=40.0)
    atm2 =  p3_atm.AtmosphereVgradient(w0=0, w1=5, h0=20, h1=60)
    figure = plt.figure(figsize=(10.24, 3.12))
    axes = figure.subplots(1, 2)
    p3_plu.plot_slice_wind_nu(atm1, n0=-100, n1=80, dn=10., e0=0., h0=-10, h1=100, dh=5., zdir=-1.,
                              show_quiver=True, show_color_bar=False,
                              title='Shear',
                              figure=figure, ax=axes[0], use_wx=True)
    p3_plu.plot_slice_wind_nu(atm2, n0=-100, n1=80, dn=10., e0=0., h0=-10, h1=100, dh=5., zdir=-1.,
                              show_quiver=True, show_color_bar=False,
                              title='Vgrad',
                              figure=figure, ax=axes[1], use_wx=True)
    plt.savefig('/tmp/atm_shear_vs_vgrad.png', dpi=120, bbox_inches='tight')
    
   
if __name__ == "__main__":
    plot_atm()
    main()
    plt.show()
