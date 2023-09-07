#! /usr/bin/env python3
'''
  Testing space indexed trajectory
'''
import math, numpy as np, matplotlib.pyplot as plt
import pdb

import pat3.vehicles.rotorcraft.multirotor_trajectory as trj
import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
import pat3.vehicles.rotorcraft.multirotor_control as ctl

class SpaceCircle:
    def __init__(self, r=1., c = [0,0], alpha0=0, dalpha=2*np.pi, ztraj=None, psitraj=None):
        self._ylen, self._nder = 4, 5
        self.r, self.c, self.alpha0, self.dalpha = r, c, alpha0, dalpha
        self.ztraj = ztraj if ztraj is not None else trj.CstOne()
        self.psitraj = psitraj if psitraj is not None else trj.CstOne()

    def get(self, l):
        Yl = np.zeros((self._ylen, self._nder))
        alpha = self.alpha0 + self.dalpha*l
        rca, rsa = self.r*np.cos(alpha), self.r*np.sin(alpha)
        # x,y
        Yl[:2,0] = self.c+np.array([ rca,  rsa])
        Yl[:2,1] = self.dalpha    *np.array([-rsa,  rca])
        Yl[:2,2] = self.dalpha**2*np.array([-rca, -rsa])
        Yl[:2,3] = self.dalpha**3*np.array([ rsa, -rca])
        Yl[:2,4] = self.dalpha**4*np.array([ rca,  rsa])
        # z, psi
        Yl[trj._z] = self.ztraj.get(l)
        Yl[trj._psi] = self.psitraj.get(l)
        return Yl

        


class SpaceIndexedTraj:
    def __init__(self, space_traj, dynamic):
        self.duration = dynamic.duration
        self._ylen, self._nder = 4, 5
        self._g = space_traj   # geometry
        self._lamdba = dynamic # dynamic

    def get(self, t):
        Yt = np.zeros((self._g._ylen, self._g._nder))
        _lambda = self._lamdba.get(t) # lambda(t), lambdadot(t)...
        _g = self._g.get(_lambda[0])  # g(lambda), dg/dlambda(lambda)...
        Yt[:,0] = _g[:,0]
        Yt[:,1] = _lambda[1]*_g[:,1]
        Yt[:,2] = _lambda[2]*_g[:,1] + _lambda[1]**2*_g[:,2]
        Yt[:,3] = _lambda[3]*_g[:,1] + 3*_lambda[1]*_lambda[2]*_g[:,2] + _lambda[1]**3*_g[:,3]
        Yt[:,4] = _lambda[4]*_g[:,1] + (3*_lambda[2]**2+4*_lambda[1]*_lambda[3])*_g[:,2] + 6*_lambda[1]**2*_lambda[2]*_g[:,3] + _lambda[1]**4*_g[:,4]
        return Yt



    
def main(duration=3., dt=1./200):
    #_trj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[1, 0, 0, 0], duration=2.)
    #_trj = trj.SmoothLine(Y00=[0, 0, 0, 0], Y10=[0.5, -0.1, 1., np.pi], duration=2.)
    #_trj = trj.Circle(v=2.)
    _trj = trj.SmoothBackAndForth(x0=[0, 0, 0, np.pi/2], x1=[1, 0, 0, np.pi/2])
    if 0:
        straj = SpaceCircle(r=1.5, c=[0,1.], alpha0=0, dalpha=2*np.pi)
        #straj = SpaceCircle(ztraj=trj.SinOne(om=12), psitraj=trj.SinOne(om=12))
        #dtraj = trj.SinOne(om=np.pi/2/duration)
        dtraj = trj.SinOne(om=np.pi/2/duration, duration=duration)
        #dtraj = trj.AffineOne(1./duration,0., duration)
        _trj = SpaceIndexedTraj(straj, dtraj)
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
    fig, axes = fdm.plot(time, Xr, window_title="State Trajectory")#, Ur)
      
    plt.show()
 

if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
