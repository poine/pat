#! /usr/bin/env python3
'''
  Testing space indexed trajectory
'''
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.vehicles.rotorcraft.multirotor_trajectory as p3_trj
import pat3.vehicles.rotorcraft.multirotor_trajectory_dev as p3_trj_dev
import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as p3_trj_fact

import pat3.vehicles.rotorcraft.multirotor_fdm as p3_fdm
import pat3.vehicles.rotorcraft.multirotor_control as p3_ctl

import pdb

def plot_dynamics(_dyn, time, fig=None, axs=None, window_title="Trajectory Dynamics"):
    lambdas = np.array([_dyn.get(t) for t in time])
    if fig is None: fig, axs = plt.subplots(p3_trj._nder, 1)
    if window_title is not None: fig.canvas.manager.set_window_title(window_title)
    fig.suptitle('$\lambda(t)$')
    for i in range(p3_trj._nder):
        axs[i].plot(time, lambdas[:,i])
    return fig, axs

def test_and_plot(traj, time, fdm, _f1, _a1, _f2, _a2, _f3, _a3):
    _f1, _a1 = plot_dynamics(traj._lamdba, time, _f1, _a1)
    Yc = np.zeros((len(time), p3_trj._ylen, p3_trj._nder))
    for i in range(0, len(time)):
        Yc[i] = traj.get(time[i])
    _f2, _a2 = p3_trj.plot(time, Yc, figure=_f2, axes=_a2)
    Xr, Ur = np.zeros((len(time), p3_fdm.sv_size)), np.zeros((len(time), p3_fdm.iv_size))
    for i in range(0, len(time)):
        Xr[i], Ur[i], Xd = p3_ctl.DiffFlatness.state_and_cmd_of_flat_output(None, Yc[i], fdm.P)
    _f3, _a3 = p3_fdm.plot(time, Xr, window_title="State Trajectory", U=None, figure=_f3, axes=_a3)
    return _f1, _a1, _f2, _a2, _f3, _a3

    
def test_duration(dt=1./200):
    _fdm = p3_fdm.MR_FDM()
    straj = p3_trj_dev.SpaceCircle(r=1.5, c=[0,1.5], alpha0=0, dalpha=2*np.pi)
    _f1, _a1, _f2, _a2, _f3, _a3 = [None]*6
    for duration in [4, 6, 8, 10]:
        time = np.arange(0, duration, dt)
        #dtraj = trj.AffineOne(1./duration,0., duration)
        dtraj = p3_trj.PolynomialOne([0,0,0,0,0], [1,0,0,0,0], duration)
        #dtraj = trj.SigmoidOne(duration, 0.5)
        _f1, _a1 = plot_dynamics(dtraj, time, _f1, _a1)
        _trj = p3_trj_dev.SpaceIndexedTraj(straj, dtraj)
        Yc = np.zeros((len(time), p3_trj._ylen, p3_trj._nder))
        for i in range(0, len(time)):
            Yc[i] = _trj.get(time[i])
        _f2, _a2 = p3_trj.plot(time, Yc, figure=_f2, axes=_a2)
        Xr, Ur = np.zeros((len(time), p3_fdm.sv_size)), np.zeros((len(time), p3_fdm.iv_size))
        for i in range(0, len(time)):
            Xr[i], Ur[i], Xd = p3_ctl.DiffFlatness.state_and_cmd_of_flat_output(None, Yc[i], _fdm.P)
        _f3, _a3 = p3_fdm.plot(time, Xr, window_title="State Trajectory", U=None, figure=_f3, axes=_a3)

        _f1.savefig('/tmp/variable_duration_lambda.png', dpi=120, bbox_inches='tight')
        _f2.savefig('/tmp/variable_duration_flat_output.png', dpi=120, bbox_inches='tight')
        _f3.savefig('/tmp/variable_duration_state.png', dpi=120, bbox_inches='tight')

import scipy.interpolate as interpolate
class SplineDyn:
    def __init__(self, duration, nb_knot=11, k=5):
        self.duration=duration
        ts = np.linspace(0, duration, nb_knot)
        lambdas = np.linspace(0, 1, nb_knot)
        #lambdas = np.array([0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.85, 0.9, 0.95, 1.])
        #lambdas = np.array([0.0000, 0.0009, 0.0196, 0.0988, 0.2666, 0.5000, 0.7334, 0.9012, 0.9804, 0.9991, 1.0000])
        self.spline = interpolate.InterpolatedUnivariateSpline(ts, lambdas, k=k)
        #bc_type=([(1, 0.), (2, 0.)], [(1, 0.)])
        #self.spline = interpolate.make_interp_spline(ts, lambdas, k=4, bc_type=bc_type)
        
    def get(self, t):
        #pdb.set_trace()
        return self.spline.derivatives(t)
        
def test_opti_dyn(duration=7., dt=1./200):
    _fdm = p3_fdm.MR_FDM()
    time = np.arange(0, duration, dt)
    traj1 = p3_trj_fact.Traj43(duration)
    d2 = SplineDyn(duration)
    traj2 = p3_trj_fact.Traj43(duration, d2)
    fig_stuff = (_f1, _a1, _f2, _a2, _f3, _a3) = [None]*6
    #for t in np.linspace(0, duration, 11):
    #    print(f'{traj1._lamdba.get(t)[0]:.4f}')

        
    fig_stuff = test_and_plot(traj1, time, _fdm, *fig_stuff)
    fig_stuff = test_and_plot(traj2, time, _fdm, *fig_stuff)

    
  
def main(duration=10., dt=1./200):
    #test_duration(dt=1./200)
    test_opti_dyn()
    plt.show()
    

if __name__ == "__main__":
    main()
