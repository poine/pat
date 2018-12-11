#! /usr/bin/env python
import sys, argparse, logging, numpy as np, matplotlib.pyplot as plt

import pdb

import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl


def parse_command_line():
    parser = argparse.ArgumentParser(description='Plot a trajectory.')
    parser.add_argument('--traj', help='the name of the trajectory', default=None)
    parser.add_argument('--list_traj', help='list all known trajectories', action='store_true', default=False)
    return parser.parse_args()

def main():
    args = parse_command_line()
    #pdb.set_trace()
    #traj = pmt.Line([0, 0, 0], [10, 10, -1], 0.1)
    if 0:
        r, v = 1., 4.; om = v/r
        traj = pmt.Circle(c=[0, 0, -1], r=1., v=4., alpha0=0., dalpha=2*np.pi, zt=pmt.SinOne(om=om))
    #traj = pmt.Circle(c=[0, 0, 1], r=1., v=4., alpha0=0., dalpha=2*np.pi, psit=pmt.CstOne(0))
    #traj = pmt.Circle(c=[0, 0, 0], r=-1., v=4., alpha0=0., dalpha=2*np.pi)
    #traj = pmt.Oval(l=2, r=1, v=4)
    #traj = pmt.DoubleOval(l=2, r=1, v=4)
    traj =  pmt.FigureOfEight()


    t0, t1, dt = 0., 2.*traj.duration, 0.01
    time = np.arange(t0, t1, dt)
    Yc = np.array([traj.get(t) for t in time])

    if 1:
        fdm = fdm.FDM()
        df = ctl.DiffFlatness()
        Xc, Uc = [], []
        for Yci in Yc:
            Xci, Uci = df.state_and_cmd_of_flat_output(Yci, fdm.P)
            Xc.append(Xci); Uc.append(Uci)
        Xc = np.array(Xc)
        Uc = np.array(Uc)
        fdm.plot(time, Xc, Uc)
        #pdb.set_trace()

    pmt.plot(time, Yc)
    pmt.plot3d(time, Yc)
    plt.show()



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()
