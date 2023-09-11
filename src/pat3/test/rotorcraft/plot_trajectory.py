#! /usr/bin/env python3
import sys, argparse, logging, numpy as np, matplotlib.pyplot as plt

import pdb

import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt
import pat3.vehicles.rotorcraft.multirotor_trajectory_factory as pmtf
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl


def parse_command_line():
    parser = argparse.ArgumentParser(description='Plot a trajectory.')
    parser.add_argument('--traj', help='the name of the trajectory', default=None)
    parser.add_argument('--repeat', help='how many times to repeat the trajectory', default=1.)
    parser.add_argument('--list', help='list all known trajectories', action='store_true', default=False)
    args = parser.parse_args()
    args.repeat = float(args.repeat)
    return args


def plot_dflat_state_and_input(time, Yc):
    _fdm = fdm.FDM()
    df = ctl.DiffFlatness()
    Xc, Uc = [], []
    for Yci in Yc:
        Xci, Uci, Xcdi = df.state_and_cmd_of_flat_output(Yci, _fdm.P)
        Xc.append(Xci); Uc.append(Uci)
    Xc = np.array(Xc)
    Uc = np.array(Uc)
    _fdm.plot(time, Xc, Uc)

def main():
    args = parse_command_line()
    if args.list:
        pmtf.print_available()
        return
    try:
        print('loading trajectory: {}'.format(args.traj))
        traj, desc = pmtf.get(args.traj)
        print('  description: {}'.format(desc))
    except KeyError: 
        print('unknown trajectory {}'.format(args.traj))
        return

    t0, t1, dt = 0., args.repeat*traj.duration, 0.01
    time = np.arange(t0, t1, dt)
    Yc = np.array([traj.get(t) for t in time])

    if 1:
        plot_dflat_state_and_input(time, Yc)

    pmt.plot(time, Yc)
    pmt.plot2d(time, Yc)
    plt.show()
   

if __name__ == '__main__':
    np.set_printoptions(linewidth=500)
    logging.basicConfig(level=logging.INFO)
    main()
