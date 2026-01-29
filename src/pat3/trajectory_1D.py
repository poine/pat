#
# 1D trajectories
#
import numpy as np
from abc import ABC, abstractmethod

# we consider 4 time derivatives, this should be a parameter...
_nder = 5

class Trajectory(ABC):
    nder = 5
    @abstractmethod
    def get(self, t): pass
    
# Constant
class CstOne(Trajectory):
    def __init__(self, c=0.):
        self.c = c
    def get(self, t):
        return np.array([self.c, 0, 0, 0, 0])

# Affine
class AffOne(Trajectory):
    def __init__(self, P0, P1):
        self.P0, self.P1 = P0, P1
        self.duration = P1[0]-P0[0]
        self.c = (P1[1]-P0[1])/self.duration

    def get(self, t):
        y = self.P0[1] + (t-self.P0[0])*self.c
        return [y, self.c, 0, 0, 0]

class AffineOne:   # very sad... merge with above
    def __init__(self, c1=-1., c2=0, duration=1.):
        self.c1, self.c2, self.duration = c1, c2, duration
        
    def get(self, t):
        return np.array([self.c1*t+self.c2, self.c1, 0, 0, 0])
 

    
# Sine
class SinOne(Trajectory):
    def __init__(self, c=0., a=1., om=1., duration=2*np.pi):
        self.duration = duration
        self.c, self.a, self.om = c, a, om
        self.t0 = 0.
        
    def get(self, t):
        alpha = self.om*(t-self.t0)
        asa, aca = self.a*np.sin(alpha), self.a*np.cos(alpha)
        return np.array([  self.c + asa,
                           self.om*aca,
                          -self.om**2*asa,
                          -self.om**3*aca,
                           self.om**4*asa   ])


# MinSnap polynomials
def arr(k,n): # arangements a(k,n) = n!/k!
    a,i = 1,n
    while i>n-k:
        a *= i; i -= 1
    return a

class PolynomialOne(Trajectory):
    def __init__(self, Y0, Y1, duration):
        self.Y0, self.Y1 = Y0, Y1
        self.duration = duration
        _der = len(Y0)    # number of time derivatives
        _order = 2*_der   # we need twice as many coefficients
        self._der, self._order = _der, _order
        # compute polynomial coefficients for time derivative zeros
        self.coefs = np.zeros((_der, _order))
        M1 = np.zeros((_der, _der))
        for i in range(_der):
            M1[i,i] = arr(i,i)
        self.coefs[0, 0:_der] = np.dot(np.linalg.inv(M1), Y0)
        M3 = np.zeros((_der, _der))
        for i in range(_der):
            for j in range(i, _der):
                M3[i,j] = arr(i,j) * duration**(j-i)
        M4 = np.zeros((_der, _der))
        for i in range(_der):
            for j in range(_der):
                M4[i,j] = arr(i, j+_der) * duration**(j-i+_der)
        M3a0k = np.dot(M3, self.coefs[0, 0:_der])
        self.coefs[0, _der:_order] = np.dot(np.linalg.inv(M4), Y1 - M3a0k)
        # fill in coefficients for the subsequent time derivatives  
        for d in range(1,_der):
            for pow in range(0,2*_der-d):
                self.coefs[d, pow] = arr(d, pow+d)*self.coefs[0, pow+d]
                
    def get(self, t):
        # Horner method for computing polynomial value
        Y = np.zeros(self._der)
        for d in range(0, self._der):
            v = self.coefs[d,-1] 
            for j in range(self._order-2, -1, -1):
                v *= t
                v += self.coefs[d,j]
                Y[d] = v
        return Y


    
#
# Simple composite (no handling of time offset)
#
class CompositeOne(Trajectory):
    def __init__(self, trajs):
        self.trajs = trajs
        self.durations = [trj.duration for trj in self.trajs]
        self.ends = np.cumsum(self.durations)
        self.duration = self.ends[-1]
        
    def get(self, t):
        cur_step = np.argmax(self.ends >= t)
        return self.trajs[cur_step].get(t)

#
# Discontinuous composite trajectory smoothed with (constant arbitrary duration) minsnap polynomials at each corner
#
class SmoothedCompositeOne(CompositeOne):
    def __init__(self, trajs, eps=0.15):
        CompositeOne.__init__(self, trajs)
        self.eps = eps
        self.corners = []
        for i in range(len(self.trajs)-1):
            Y0, Y1 = self.trajs[i].get(self.ends[i]-self.eps), self.trajs[i+1].get(self.ends[i]+self.eps)
            self.corners.append(PolynomialOne(Y0, Y1, 2*self.eps))
        
    def get(self, t):
        cur_step = np.argmax(self.ends >= t)
        if cur_step<len(self.trajs)-1 and t>self.ends[cur_step]-self.eps:
            return self.corners[cur_step].get(t-self.ends[cur_step]+self.eps)
        elif cur_step > 0 and t<self.ends[cur_step-1]+self.eps:
            return self.corners[cur_step-1].get(t-self.ends[cur_step-1]+self.eps)
        else: return self.trajs[cur_step].get(t)

#
# Simple constant velocity spatial index with zero velocity at start and stop
# (used for optimization)
#
class SmoothStopStopCstVel(Trajectory):
    def __init__(self, dt_acc, dl_acc, dt_cruise):
        self.dt_acc, self.dl_acc, self.dt_cruise = dt_acc, dl_acc, dt_cruise
        self.duration = 2*dt_acc+dt_cruise
        self.c = (1.-2*dl_acc)/dt_cruise # cruise slope
        self.intro = PolynomialOne([0, 0, 0, 0, 0], [dl_acc, self.c, 0, 0, 0], dt_acc)
        self.outro = PolynomialOne([1-dl_acc, self.c, 0, 0, 0], [1, 0, 0, 0, 0], dt_acc)
       
    def get(self, t):
        if t < self.dt_acc: return self.intro.get(t)
        elif t < self.dt_acc+self.dt_cruise: return [self.dl_acc+(t-self.dt_acc)*self.c, self.c, 0, 0, 0]
        else: return self.outro.get(t-self.dt_acc-self.dt_cruise)
