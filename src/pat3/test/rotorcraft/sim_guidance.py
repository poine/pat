#! /usr/bin/env python3

import sys, argparse, math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.algebra as pal, pat3.utils as pmu, pat3.atmosphere as patm
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl
import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt
import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as pmtf

def parse_command_line():
    parser = argparse.ArgumentParser(description='Plot a trajectory.')
    parser.add_argument('--traj', help='the name of the trajectory', default=None)
    parser.add_argument('--list', help='list all known trajectories', action='store_true', default=False)
    return parser.parse_args()

def main(dt=0.005):
    args = parse_command_line()
    if args.list:
        print('available trajectories:')
        for n in pmtf.list():
            print(' {}'.format(n))
    try:
        print('loading trajectory: {}'.format(args.traj))
        _traj, _desc = pmtf.get(args.traj)
        print('  desc: {}'.format(_desc))
    except KeyError: 
        print('unknown trajectory {}'.format(args.traj))
        return
    _fdm = fdm.MR_FDM()
    #_fdm = fdm.UFOFDM()
    #_fdm = fdm.SolidFDM()
    _ctl_input = ctl.TrajRef(_traj, _fdm)
    _ctl = ctl.PosController(_fdm, _ctl_input)
    _atm = np.array([0,0,0])#patm.AtmosphereCstWind()#FIXME
    sim = pmu.Sim(_fdm, _ctl, _atm)

    time = np.arange(0, _traj.duration, dt)
    X, U = np.zeros((len(time), fdm.sv_size)), np.zeros((len(time), _fdm.iv_size))
    Xref = np.zeros((len(time), _ctl.ref_size))
    X[0] = sim.reset(time[0], _ctl_input.get(time[0])[2])
    for i in range(1, len(time)):
        U[i-1], X[i] = sim.run(time[i])
        Xref[i-1] = _ctl.Xref
    U[-1], Xref[-1] = U[-2], Xref[-2]
        
    figure = _fdm.plot(time, X, U)
    _ctl.plot(time, U, Xref)
    plt.show()
    
if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    #_traj =  pmt.Circle(v=1.)
    #_traj =  pmt.FigureOfEight(v=1.)
    #_traj, _desc = pmtf.get('circle_with_intro_slow')
    main()
