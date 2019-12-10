#! /usr/bin/env python
import os, sys, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import control.matlab

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil
import pat3.utils as p3_u
import pat3.atmosphere as p3_atm
import pat3.frames as p3_fr
import pat3.plot_utils as p3_pu

import test_02_att_ctl

#
# I am runnin a few simulations in order to find the best circle radius
# for climbing in a thermal
#
def find_best_radius(dm, atm, va=9., compute=False, plot_atm=False, plot_traj=False):
    if plot_atm:
        p3_pu.plot_slice_wind(atm, xmax=100, dx=2.5)
        plt.show()

    savefile_name = '/tmp/thermal_optim_{:.1f}.npz'.format(va)
    if compute:
        time, Xe, Ue, phi_sp, theta_sp = test_02_att_ctl.get_sim_defaults(dm, tf=15, trim_args = {'h':0, 'va':va, 'throttle':0})
        vzs, phis = [], []
        for _phi_sp in np.deg2rad(np.arange(0, 45, 1.)):
            phi_sp = _phi_sp*np.ones(len(time))
            time, X, U = test_02_att_ctl.run_simulation(dm, time, Xe, Ue, phi_sp, theta_sp, plot=plot_traj)
            if plot_traj: plt.show()
            Xee = dm.state_six_dof_euclidian_euler(X[-1], atm=None)
            vzs.append(Xee[p3_fr.SixDOFEuclidianEuler.sv_zd])
            phis.append(X[-1][dm.sv_phi])
            #print('vz {:.2f}m/s'.format(Xee[p3_fr.SixDOFEuclidianEuler.sv_zd]))
            #print('phi {}'.format(np.rad2deg(X[-1][dm.sv_phi])))
        phis, vzs = np.array(phis), np.array(vzs)
        g=9.81 
        Rs = va**2/g/np.tan(phis)
        #atm = p3_atm.AtmosphereThermal1()
        updrafts = np.array([atm.get_wind(pos_ned=[_r, 0, 0], t=0) for _r in Rs])
        np.savez(savefile_name, vzs=vzs, phis=phis, Rs=Rs, updrafts=updrafts)
    else:
        _data =  np.load(savefile_name)
        vzs, phis, Rs, updrafts = [_data[k] for k in ['vzs', 'phis', 'Rs', 'updrafts']]
        
    climb_rate = updrafts[1:,2]+vzs[1:]
    best_climb_idx = np.argmin(climb_rate)
    res_txt = 'best vz {:.2f} radius {:.1f} roll {:.1f}'.format(climb_rate[best_climb_idx], Rs[best_climb_idx], np.rad2deg(phis[best_climb_idx]))
    print('va {:.1f} '.format(va)+ res_txt)
    good_Rs = Rs < 100.
    fig = p3_pu.prepare_fig(window_title='va = {:.1f} m/s'.format(va))
    ax = plt.subplot(1, 3, 1)
    plt.plot(np.rad2deg(phis), vzs, label='aircraft')
    plt.plot(np.rad2deg(phis), updrafts[:,2], label='updraft')
    p3_pu.decorate(ax, title='Vz', xlab='roll in deg', ylab='vz in m/s', legend=True)
    plt.subplot(1, 3, 2)
    if 1:
        plt.plot(np.rad2deg(phis[good_Rs]), Rs[good_Rs])
        p3_pu.decorate(plt.gca(), title='Radius', xlab='roll in deg', ylab='radius in m')
    if 0:
        plt.plot(np.rad2deg(phis), 1/Rs)
        p3_pu.decorate(plt.gca(), title='Curvature', xlab='roll in deg', ylab='curvature in m^-1')
    #pdb.set_trace()
    ax = plt.subplot(1, 3, 3)
    #plt.plot(Rs[good_Rs], vzs[good_Rs], label='aircraft')
    #plt.plot(Rs[good_Rs], updrafts[good_Rs,2], label='updrafts')
    plt.plot(Rs[1:], vzs[1:], label='aircraft')
    plt.plot(Rs[1:], updrafts[1:,2], label='updrafts')
    plt.plot(Rs[1:], updrafts[1:,2]+vzs[1:], label='sum')
    p3_pu.decorate(ax, title='Vz\n{}'.format(res_txt), xlab='R in m', ylab='vz in m/s', legend=True)



##
## computing FFTs on Vz
##
def pts_on_circle(alpha, _c, _r): return [_c[0]+_r*np.cos(alpha), _c[1]+_r*np.sin(alpha)]

def _foo(alphas, wc1, dt, time, _c):
    w_mean, w_ampl = np.mean(wc1), np.max(wc1)-np.min(wc1)
    spectrum = np.fft.fft(wc1-w_mean)
    #spectrum = np.fft.fft(wc1)
    freqs = np.fft.fftfreq(spectrum.shape[-1], dt)
    magnitude = np.abs(spectrum)
    max_mag_idx =  np.argmax(magnitude)
    #threshold = np.max(magnitude)/10000
    #spectrum2 = np.array(spectrum)
    #spectrum2[magnitude<threshold] = 0.
    phases = np.angle(spectrum)
    
    #plt.plot(1./freq, spectrum.real, 1./freq, spectrum.imag)
    #plt.plot(1./freq, phase)
    periods = 1./freqs 
    if 1:
        ax=plt.subplot(2,3,3)
        plt.plot(periods, magnitude)
        period = periods[max_mag_idx]
        print('period {:.1f}s'.format(period))
    if 1:
        ax=plt.subplot(2,3,6)
        phase = phases[max_mag_idx]
        phases[magnitude < 200] = 0
        plt.plot(periods, phases)
        print('phase {:.1f} deg'.format(np.rad2deg(phase)))
        #plt.subplot(2,3,5)
        #plt.plot(time, w_mean+w_ampl/2*np.sin(time/period*2*np.pi+phase+np.pi/2))
        ax=plt.subplot(1,3,1)
        p1, p2 = _c, pts_on_circle(np.pi/2-phase, _c, 15)
        #plt.plot([p1[0], p2[0]], [p1[1], p2[1]])

def plot_gradient(atm, xmax=50, dx=5.):
    h=0
    xlist, ylist = np.arange(-xmax, xmax, dx), np.arange(-xmax, xmax, dx)
    x, y = np.meshgrid(xlist, ylist)
    wz = np.zeros_like(y)
    nx, ny = x.shape
    for ix in range(x.shape[0]):
        for iy in range(x.shape[1]):
            wz[ix, iy] = -atm.get_wind([x[ix, iy], y[ix, iy], -h], t=0)[2]

    #_c, _r = [10, 0], 20.
    #_cs = [[-10, 0], [0, 0], [10, 0], [20, 0], [30, 0]]
    _cs = [[-10, 0], [10, 0], [30, 0]]
    #_cs = [[-40, 0]]
    ax=plt.subplot(1,3,1)
    ax.axis('equal')
    cp = ax.contourf(x, y, wz, alpha=0.5)
    plt.gcf().colorbar(cp)
    ax.set_title('Updraft')
    ax=plt.subplot(2,3,2)
    p3_pu.decorate(ax, title='Vz', xlab='alpha in deg', ylab='vz in m/s', legend=True)
    ax=plt.subplot(2,3,5)
    p3_pu.decorate(ax, title='Vz', xlab='time in s', ylab='vz in m/s', legend=True)
    for _c in _cs:
        _r = 20.
        va, dt = 9., 0.01
        time = np.arange(0, 40., dt)
        alphas = va*time/_r#np.arange(0, 4*np.pi, 0.1)
        c_pts = np.array([pts_on_circle(_alpha, _c, _r) for _alpha in alphas])
        ax=plt.subplot(1,3,1)
        plt.plot(c_pts[:,0], c_pts[:,1])
        ax=plt.subplot(2,3,2)
        wc1 = [-atm.get_wind([_x, _y, -h], t=0)[2] for _x, _y in zip(c_pts[:,0], c_pts[:,1])]
        plt.plot(np.rad2deg(alphas), wc1)
        ax=plt.subplot(2,3,5)
        plt.plot(time, wc1)
        ax=plt.subplot(2,3,6)
        _foo(alphas, wc1, dt, time, _c)

def test_circles(atm, xmax=50, dx=5.):
    pass

        
def main(param_filename):
    dm = p1_fw_dyn.DynamicModel(param_filename)
    atm = p3_atm.AtmosphereThermal1()
    #atm.set_params(xc=0, yc=10, zi=1500, wstar=300)
    #find_best_radius(dm, atm, va=9., compute=True, plot_atm=False, plot_traj=True)
    plot_gradient(atm)
    plt.show()  
    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
    main(param_filename)

