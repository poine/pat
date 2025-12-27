import numpy as np

import pat3.vehicles.rotorcraft.multirotor_trajectory as trj
#import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
#import pat3.vehicles.rotorcraft.multirotor_control as ctl


class SpaceIndexedLine(trj.Line):
    def __init__(self, p1, p2, psi):
        self._ylen, self._nder = 4, 5
        trj.Line.__init__(self, p1, p2, 1/np.linalg.norm(p2-p1), psi)

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

import scipy.interpolate as interpolate
class SpaceWaypoints:
    def __init__(self, waypoints):
        self._ylen, self._nder = 4, 5
        self.waypoints = np.array(waypoints)
        l = np.linspace(0, 1, len(self.waypoints))
        self.splines = [interpolate.InterpolatedUnivariateSpline(l, self.waypoints[:,i], k=4) for i in range(3)]
        
    def get(self, l):
        Yl = np.zeros((self._ylen, self._nder))
        Yl[0:3] = [self.splines[i].derivatives(l) for i in range(3)]
        return Yl
 
class SpaceIndexedTraj:
    def __init__(self, geometry, dynamic):
        self.duration = dynamic.duration
        self._ylen, self._nder = 4, 5
        self._geom, self._dyn = geometry, dynamic
        self.t0 = 0.

    def set_dyn(self, dyn): self._dyn = dyn

    def get(self, t):
        Yt = np.zeros((self._geom._ylen, self._geom._nder))
        _lambda = self._dyn.get(t) # lambda(t), lambdadot(t)...
        _lambda[0] = np.clip(_lambda[0], 0., 1.)   # protect ourselvf against unruly dynamics 
        _g = self._geom.get(_lambda[0])  # g(lambda), dg/dlambda(lambda)...
        Yt[:,0] = _g[:,0]
        Yt[:,1] = _lambda[1]*_g[:,1]
        Yt[:,2] = _lambda[2]*_g[:,1] + _lambda[1]**2*_g[:,2]
        Yt[:,3] = _lambda[3]*_g[:,1] + 3*_lambda[1]*_lambda[2]*_g[:,2] + _lambda[1]**3*_g[:,3]
        Yt[:,4] = _lambda[4]*_g[:,1] + (3*_lambda[2]**2+4*_lambda[1]*_lambda[3])*_g[:,2] + 6*_lambda[1]**2*_lambda[2]*_g[:,3] + _lambda[1]**4*_g[:,4]
        return Yt




