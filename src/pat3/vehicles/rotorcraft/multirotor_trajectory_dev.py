import numpy as np

import pat3.vehicles.rotorcraft.multirotor_trajectory as trj
#import pat3.vehicles.rotorcraft.multirotor_fdm as fdm
#import pat3.vehicles.rotorcraft.multirotor_control as ctl

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
        print(l)
        
    def get(self, l):
        Yl = np.zeros((self._ylen, self._nder))
        Yl[0:3] = [self.splines[i].derivatives(l) for i in range(3)]
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




class FooTraj(SpaceIndexedTraj):
    def __init__(self, duration=6):
        r, v = 1.5, 2.; om = v/r; alpha0 = 0
        psit = trj.AffineOne(om, alpha0)
        straj = SpaceCircle(r=r, c=[0,1.], alpha0=0, dalpha=4*np.pi, psitraj=psit)
        #dtraj = trj.SinOne(om=np.pi/2/duration, duration=duration)
        dtraj = trj.PolynomialOne([0,0,0,0,0], [1, 0, 0, 0, 0], duration=10.)
        #dtraj = trj.AffineOne(1./duration,0., duration)
        
        SpaceIndexedTraj.__init__(self, straj, dtraj)
