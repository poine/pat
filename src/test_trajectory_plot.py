#! /usr/bin/env python
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
    return parser.parse_args()

def main():
    args = parse_command_line()
    if args.list:
        print('available trajectories:')
        for n in pmtf.list():
            print(' {}'.format(n))
        return

    try:
        print('loading trajectory: {}'.format(args.traj))
        traj, desc = pmtf.get(args.traj)
        print('  description: {}'.format(desc))
        
        t0, t1, dt = 0., args.repeat*traj.duration, 0.01
        time = np.arange(t0, t1, dt)
        Yc = np.array([traj.get(t) for t in time])

        if 1:
            _fdm = fdm.FDM()
            df = ctl.DiffFlatness()
            Xc, Uc = [], []
            for Yci in Yc:
                Xci, Uci = df.state_and_cmd_of_flat_output(Yci, _fdm.P)
                Xc.append(Xci); Uc.append(Uci)
            Xc = np.array(Xc)
            Uc = np.array(Uc)
            _fdm.plot(time, Xc, Uc)
            #pdb.set_trace()

        pmt.plot(time, Yc)
        pmt.plot3d(time, Yc)
        plt.show()
    except KeyError: 
        print('unknown trajectory {}'.format(args.traj))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
