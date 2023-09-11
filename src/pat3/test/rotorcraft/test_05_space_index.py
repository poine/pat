#! /usr/bin/env python3
'''
  Testing space indexed trajectory
'''
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.vehicles.rotorcraft.multirotor_trajectory as trj
import pat3.vehicles.rotorcraft.multirotor_trajectory_dev as trj_dev


import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl

import pdb

class WaypointTraj:
    def __init__(self):
        self.duration = 1.
        self.waypoints = np.array([[0,0,0,0], [1, 0, 0, 0], [1, 1, 0, 0], [2, 1, 0, 0]])
        self.nb_segments = self.waypoints.shape[0] - 1
        self.segment_duration = self.duration/self.nb_segments
        self.segments = []
        for i in range(self.nb_segments):
            x0 = [self.waypoints[i,0], 0, 0, 0, 0]
            x1 = [self.waypoints[i+1,0], 0, 0, 0, 0]
            y0 = [self.waypoints[i,1], 0, 0, 0, 0]
            y1 = [self.waypoints[i+1,1], 0, 0, 0, 0]
            _px = trj.PolynomialOne(x0, x1, self.segment_duration)
            _py = trj.PolynomialOne(y0, y1, self.segment_duration)
            self.segments.append([_px, _py])
        
    def reset(self, t0): self.t0 = t0
        
    def get(self, t):
        Yc = np.zeros((4,5))
        id_segment = int(np.floor(t*self.nb_segments))
        _t = t - id_segment*self.segment_duration
        Yc[0] = self.segments[id_segment][0].get(_t)
        Yc[1] = self.segments[id_segment][1].get(_t)
        print(t, _t, id_segment, Yc)
        
        return Yc

    
def main(duration=10., dt=1./200):
    if 1:
        #straj = trj_dev.SpaceCircle(r=1.5, c=[0,1.], alpha0=0, dalpha=2*np.pi)
        straj = trj_dev.SpaceWaypoints([[0, 0, 0], [1, 1, 0.5], [0, 2, 0], [1, 3, 0.5], [0, 4, 0]])
        #straj = SpaceCircle(ztraj=trj.SinOne(om=12), psitraj=trj.SinOne(om=12))
        #dtraj = trj.SinOne(om=np.pi/2/duration)
        #dtraj = trj.SinOne(om=np.pi/2/duration, duration=duration)
        dtraj = trj.AffineOne(1./duration,0., duration)
        _trj = trj_dev.SpaceIndexedTraj(straj, dtraj)
    else:
        #_trj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[1, 0, 0, 0], duration=2.)
        #_trj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[0.5, -0.1, 1., np.pi], duration=2.)
        #_trj = trj.Circle(v=2.)
        #_trj = trj.SmoothBackAndForth(x0=[0, 0, 0, np.pi/2], x1=[1, 0, 0, np.pi/2])
        _trj = WaypointTraj()
        
    time = np.arange(0, _trj.duration, dt)
    Yc = np.zeros((len(time), trj._ylen, trj._nder))
    for i in range(0, len(time)):
        Yc[i] = _trj.get(time[i])

    trj.check_consistency(time, Yc)
    #trj.plot(time, Yc)
    trj.plot2d(time, Yc)

    _fdm = fdm.MR_FDM()
    Xr, Ur = np.zeros((len(time), fdm.sv_size)), np.zeros((len(time), fdm.iv_size))
    for i in range(0, len(time)):
        Xr[i], Ur[i], Xd = ctl.DiffFlatness.state_and_cmd_of_flat_output(None, Yc[i], _fdm.P)
    #fig, axes = fdm.plot(time, Xr, window_title="State Trajectory")#, Ur)
      
    plt.show()
 

if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
