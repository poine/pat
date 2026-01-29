#-*- coding: utf-8 -*-

import math, numpy as np, matplotlib.pyplot as plt, matplotlib

'''
  (Output) Trajectories for a quadrotor
  Y = [x, y, z, psi] plus their 4 time derivatives
  Y constitutes a flat output and can hence be used to compute the state and input
  of the quadrotor (see vehicle/rotorcraft/multirotor_control.py)
'''
import pat3.plot_utils as ppu, pat3.algebra as pal
import pat3.utils as p3_u
import pat3.trajectory_1D as p3_t1D

_x, _y, _z, _psi, _ylen = range(5)
# we consider 4 time derivatives, this should be a parameter...
_nder = 5


# ALL 1D trajectories moved to pat3.trajectory_1D 
#
# Scalar trajectories (aka 1D)
# 
# class CstOne:
#     def __init__(self, c=-1.):
#         self.c = c
#     def get(self, t):
#         return np.array([self.c, 0, 0, 0, 0])

# class AffineOne:
#     def __init__(self, c1=-1., c2=0, duration=1.):
#         self.c1, self.c2, self.duration = c1, c2, duration
        
#     def get(self, t):
#         return np.array([self.c1*t+self.c2, self.c1, 0, 0, 0])
    
# class SinOne:
#     def __init__(self, c=0., a=1., om=1., duration=2*np.pi):
#         self.duration = duration
#         self.c, self.a, self.om = c, a, om
#         self.t0 = 0.
        
#     def get(self, t):
#         alpha = self.om*(t-self.t0)
#         asa, aca = self.a*np.sin(alpha), self.a*np.cos(alpha)
#         return np.array([  self.c + asa,
#                            self.om*aca,
#                           -self.om**2*asa,
#                           -self.om**3*aca,
#                            self.om**4*asa   ])

# # TODO - compute derivatives
# class SigmoidOne:
#     def __init__(self, duration=5, l=2):
#         self.duration = duration
#         self.l = l
        
#     def get(self, t):
#         a = np.tan(-np.pi/2+t*np.pi/self.duration)
#         b = 1./(1+np.exp(-self.l*a))
#         return [b, 0, 0, 0, 0]
        

    
# def arr(k,n):
#     '''arangements a(k,n) = n!/k!'''
#     a,i = 1,n
#     while i>n-k:
#         a *= i
#         i -= 1
#     return a

# class PolynomialOne:
#     def __init__(self, Y0, Y1, duration):
#         self.Y0, self.Y1 = Y0, Y1
#         self.duration = duration
#         _der = len(Y0)    # number of time derivatives
#         _order = 2*_der   # we need twice as many coefficients
#         self._der, self._order = _der, _order
#         # compute polynomial coefficients for time derivative zeros
#         self.coefs = np.zeros((_der, _order))
#         M1 = np.zeros((_der, _der))
#         for i in range(_der):
#             M1[i,i] = arr(i,i)
#         self.coefs[0, 0:_der] = np.dot(np.linalg.inv(M1), Y0)
#         M3 = np.zeros((_der, _der))
#         for i in range(_der):
#             for j in range(i, _der):
#                 M3[i,j] = arr(i,j) * duration**(j-i)
#         M4 = np.zeros((_der, _der))
#         for i in range(_der):
#             for j in range(_der):
#                 M4[i,j] = arr(i, j+_der) * duration**(j-i+_der)
#         M3a0k = np.dot(M3, self.coefs[0, 0:_der])
#         self.coefs[0, _der:_order] = np.dot(np.linalg.inv(M4), Y1 - M3a0k)
#         # fill in coefficients for the subsequent time derivatives  
#         for d in range(1,_der):
#             for pow in range(0,2*_der-d):
#                 self.coefs[d, pow] = arr(d, pow+d)*self.coefs[0, pow+d]
                
#     def get(self, t):
#         # Horner method for computing polynomial value
#         Y = np.zeros(self._der)
#         for d in range(0, self._der):
#             v = self.coefs[d,-1] 
#             for j in range(self._order-2, -1, -1):
#                 v *= t
#                 v += self.coefs[d,j]
#                 Y[d] = v
#         return Y
        

# import scipy.interpolate as interpolate
# class SplinesOne:
#     def __init__(self, points_t, points_l):
#         self.duration = points_t[-1]
#         self._nder = 5
#         self.p_t, self.p_l = points_t, points_l
#         self.spline = interpolate.InterpolatedUnivariateSpline(points_t, points_l, k=4)
        
#     def get(self, t):
#         return self.spline.derivatives(t)




    
from abc import ABC, abstractmethod
#
# Solid Trajectories (aka 3D)
# 

class Trajectory(ABC):
    def is_space_indexed(self): return False  # space indexed trajectories return True
    def has_dyn_ctl_pts(self): return False   # space_indexed trajectories might have them  
    def has_waypoints(self): return False     # 
    @abstractmethod
    def get(self, t): pass
    
class Cst(Trajectory):
    def __init__(self, Y00, duration=1.):
        self.duration = duration
        self.Yc = np.zeros((_ylen, _nder))
        self.Yc[:,0] = Y00
        
    def reset(self, t0): pass

    def get(self, t):
        return self.Yc
        


class Circle(Trajectory):

    def __init__(self, c=[0, 0, 0], r=1., v=2., alpha0=0, dalpha=2*np.pi, zt=None, psit=None):
        self.c, self.r, self.v = np.asarray(c), r, v # center, radius, velocity
        self.alpha0, self.dalpha = alpha0, dalpha    # start angle, angle span
        self.omega = self.v/self.r
        self.length = dalpha*np.abs(r)
        self.duration = self.length/v
        self.zt = zt if zt is not None else p3_t1D.CstOne(c[_z])
        self.psit = psit if psit is not None else p3_t1D.AffineOne(self.omega, self.alpha0+np.sign(self.r)*np.pi/2)
        #self.psit = psit if psit is not None else AffineOne(self.omega, self.alpha0)
        self.t0 = 0.
        
    def reset(self, t0):
        self.t0 = t0
        #print('in circle reset {}'.format(t0))
     
    def get(self, t):
        dt = t-self.t0
        alpha = self.omega*(dt) + self.alpha0
        #print('in circle get t:{} dt: {} alpha {}'.format(t, dt, alpha))
        rca, rsa = np.abs(self.r)*np.cos(alpha), np.abs(self.r)*np.sin(alpha) 
        Yc = np.zeros((5,4))
        ## xy
        Yc[0,:_z] = self.c[:_z] + [rca, rsa]
        Yc[1,:_z] = [-self.omega   *rsa,  self.omega   *rca]
        Yc[2,:_z] = [-self.omega**2*rca, -self.omega**2*rsa]
        Yc[3,:_z] = [ self.omega**3*rsa, -self.omega**3*rca]
        Yc[4,:_z] = [ self.omega**4*rca,  self.omega**4*rsa]
        ## z
        Yc[:,_z] = self.zt.get(dt)
        ## psi
        Yc[:,_psi] =  self.psit.get(dt)
        Yc[0,_psi]  = pal.norm_mpi_pi(Yc[0,_psi])
        return Yc.T



class Line:
    def __init__(self, p1, p2, v=2., psi=None):
        self.p1, self.p2, self.v = np.asarray(p1), np.asarray(p2), v # ends and velocity 
        dep = self.p2-self.p1
        self.length = np.linalg.norm(dep)   # length
        self.un = dep/self.length           # unit vector
        self.psi = psi if psi is not None else np.arctan2(self.un[1], self.un[0])
        self.duration = self.length/self.v  # duration
        self.t0 = 0.

    def reset(self, t0): self.t0 = t0
        
    def get(self, t):
        Yc = np.zeros((_nder, _ylen))
        Yc[0,:3] = self.p1 + self.un*self.v*(t-self.t0)
        Yc[1,:3] =           self.un*self.v
        Yc[0,3] = self.psi
        return Yc.T



class Ellipse:
    pass


class SmoothLine:
    def __init__(self, Y00=[0, 0, 0, 0], Y10=[1, 0, 0, 0], duration=1.):
        self.duration = duration
        Y0 = np.zeros((_ylen, _nder))# [_x, _y, _z, _psi, _ylen]
        if len(np.asarray(Y00).shape) == 1: # we only got zero order derivatives
            Y0[:,0] = Y00
        else:
            Y0 = Y00
        Y1 = np.zeros((_ylen, _nder))
        if len(np.asarray(Y10).shape) == 1: # we only got zero order derivatives
            Y1[:,0] = Y10
        else:
            Y1 = Y10
        if Y1[_psi,0]-Y0[_psi,0] > np.pi: Y1[_psi,0]-=2*np.pi
        if Y1[_psi,0]-Y0[_psi,0] < -np.pi: Y1[_psi,0]+=2*np.pi
        
        self._polys = [p3_t1D.PolynomialOne(Y0[i], Y1[i], self.duration) for i in range(_ylen)]
        self.t0 = 0
        
    def reset(self, t0):
        self.t0 = t0
        
    def get(self, t):
        Yc = np.array([p.get(t-self.t0) for p in self._polys])
        return Yc
        
    
    
class CompositeTraj(Trajectory):
    def __init__(self, steps):
        self.steps = steps
        self.steps_dur = [s.duration for s in self.steps]
        self.steps_end = np.cumsum(self.steps_dur)
        self.duration = np.sum(self.steps_dur)
        for s, st in zip(self.steps[1:], self.steps_end):
            s.reset(st)
        self.t0 = 0.   

    def reset(self, t0):
        self.t0 = t0
        
    def get(self, t):
        dt = t - self.t0
        Yc = np.zeros((5,4))
        dt_lapse = math.fmod(dt, self.duration)
        cur_step = np.argmax(self.steps_end > dt_lapse)
        Yc = self.steps[cur_step].get(dt_lapse)
        return Yc
    
class Oval(CompositeTraj):
    def __init__(self, l, r, v=2., z=-1):
        self.l, self.r, self.v = l, r, v # length, radius and velocity 

        c1, c2 = np.array([-l, 0, z]), np.array([l, 0, z])
        p1, p2 = np.array([-l, -r, z]), np.array([l, -r, z])
        p3, p4 = np.array([l, r, z]), np.array([-l, r, z])
        steps = [Line(p1, p2, v),
                 Circle(c2, r, v, -np.pi/2, np.pi),
                 Line(p3, p4, v),
                 Circle(c1, r, v, np.pi/2, np.pi)]
        CompositeTraj.__init__(self, steps)
    

class DoubleOval(CompositeTraj):
    def __init__(self, l, r, v=2., z=-1):
        p1, p2 = np.array([0, 0, z]), np.array([l, 0, z])
        p3, p4 = np.array([l, 2*r, z]), np.array([0, 2*r, z])
        p5, p6 = np.array([0, 4*r, z]), np.array([l, 4*r, z])
        c1, c2 = np.array([l, r, z]), np.array([0, 3*r, z])
        c3, c4 = np.array([l, 3*r, z]), np.array([0, r, z])
        steps = [Line(p1, p2, v),
                 Circle(c1, r,  v, -np.pi/2, np.pi),
                 Line(p3, p4, v),
                 Circle(c2, -r, v, -np.pi/2, np.pi),
                 Line(p5, p6, v),
                 Circle(c3, -r, v,  np.pi/2, np.pi),
                 Line(p3, p4, v),
                 Circle(c4, r,  v,  np.pi/2, np.pi),
        ] 
        CompositeTraj.__init__(self, steps)
        
        
class FigureOfEight(CompositeTraj):
    def __init__(self, r=1., v=2., z=-0.25):
        c1, c2 = np.array([r, 0, z]), np.array([-r, 0, z])
        steps = [Circle(c1,  r,  v, -np.pi, 2*np.pi),
                 Circle(c2, -r,  v, 0,      2*np.pi) ]
        CompositeTraj.__init__(self, steps)
          
class SmoothBackAndForth(CompositeTraj):
    def __init__(self, Y0=[0, 0, 0.5, 0], Y1=[1, 0, -0.5, 0], dt_move=2., dt_stay=1.):
        steps = [SmoothLine(Y0, Y1, duration=dt_move),
                 Cst(Y1, dt_stay),
                 SmoothLine(Y1, Y0, duration=dt_move),
                 Cst(Y0, dt_stay)]
        CompositeTraj.__init__(self, steps)     

class CircleWithIntro(CompositeTraj):
    
    def __init__(self, Y0=[0, 0, -1.5, 0], c=[0, 0, -1.5], r=1., v=3., dt_intro=1.8, dt_stay=0., psit=None):
        eps = np.deg2rad(2)
        circle = Circle(c,  r,  v, np.pi/2-eps, 2*np.pi+2*eps, psit=psit)
        #Y0 = [0, 0, -1.5, 0]
        Y1 = circle.get(0.)
        Y2 = circle.get(circle.duration)
        steps = [SmoothLine(Y0, Y1, duration=dt_intro),
                 circle,
                 SmoothLine(Y2, Y0, duration=dt_intro),
                 Cst(Y0, dt_stay)] 
        CompositeTraj.__init__(self, steps)
        

class RefModelTraj:
    def __init__(self, traj_setpoint, dt=0.01, Y0=None):
        dyns = [[8., 0.7, 4., 0.9], [8., 0.7, 4., 0.9], [8., 0.7, 4., 0.9], [5., 0.7, 2., 0.9]]
        self.refs = [p3_u.FourthOrdLinRef(om1, xi1, om2, xi2) for om1, xi1, om2, xi2 in dyns]
        self.duration = traj_setpoint.duration
        self.time = np.arange(0, self.duration, dt)
        self.Ysp = np.array([traj_setpoint.get(_t) for _t in self.time])
        for i in range(4):
            self.refs[i].reset(self.Ysp[0,i] if Y0 is None else Y0[i])

        self.Ys = np.zeros_like(self.Ysp)
        for i in range(1, len(self.time)):
            for ycmp in range(4):
                if ycmp==3:  # wrap pesky heading angle
                    err = self.Ys[i-1, ycmp, 0] - self.Ysp[i, ycmp, 0]
                    if err > np.pi: self.Ysp[i, ycmp, 0] += 2*np.pi
                    elif err < -np.pi: self.Ysp[i, ycmp, 0] -= 2*np.pi
                self.Ys[i, ycmp] = self.refs[ycmp].run(dt, self.Ysp[i, ycmp, 0])

        
    def reset(self, t0): pass
        
    def get(self, t):
        idx = np.argmax(self.time > t)
        Yc = self.Ys[idx]
        return Yc



# check trajectory consistency
def check_consistency(time, Y):
    Ycheck = np.zeros_like(Y)
    # compute numerical differentiation of provided trajectory
    Ycheck[:,:,1] = np.gradient(Y[:,:,0], time[1]-time[0], axis=0)
    # compute further numerical differentiations
    for j in range(2, _nder):
        Ycheck[:,:,j] = np.gradient(Ycheck[:,:,j-1], time[1]-time[0], axis=0)
    
    figure, axes = plot(time, Y)
    _s = 4
    for i in range(_ylen): # x, y, z, psi
        for j in range(1, _nder): # the four time derivatives
            axes[j,i].plot(time[j:-j], np.rad2deg(Ycheck[j:-j,i,j]) if i == _psi else Ycheck[j:-j,i,j], label="check")

    
def plot(time, Yc, figure=None, axes=None, window_title="Flat Output Trajectory"):
    figure = ppu.prepare_fig(figure, window_title, (20.48, 10.24))
    plots = [("$x$",       "m",     0.5, Yc[:,_x, 0]),
             ("$y$",       "m",     0.5, Yc[:,_y, 0]),
             ("$z$",       "m",     0.5, Yc[:,_z, 0]),
             ("$\\psi$",    "deg",   0.5, np.rad2deg(Yc[:,_psi, 0])),
             ("$x^{(1)}$", "m/s",   0.5, Yc[:,_x, 1]),
             ("$y^{(1)}$", "m/s",   0.5, Yc[:,_y, 1]),
             ("$z^{(1)}$", "m/s",   0.5, Yc[:,_z, 1]),
             ("$\\psi^{(1)}$", "deg/s",   0.5, np.rad2deg(Yc[:,_psi, 1])),
             ("$x^{(2)}$", "m/s2",  0.5, Yc[:,_x, 2]),
             ("$y^{(2)}$", "m/s2",  0.5, Yc[:,_y, 2]),
             ("$z^{(2)}$", "m/s2",  0.5, Yc[:,_z, 2]),
             ("$\\psi^{(2)}$", "deg/s2",   0.5, np.rad2deg(Yc[:,_psi, 2])),
             ("$x^{(3)}$", "m/s3",  0.5, Yc[:,_x, 3]),
             ("$y^{(3)}$", "m/s3",  0.5, Yc[:,_y, 3]),
             ("$z^{(3)}$", "m/s3",  0.5, Yc[:,_z, 3]),
             ("$\\psi^{(3)}$", "deg/s3",   0.5, np.rad2deg(Yc[:,_psi, 3])),
             ("$x^{(4)}$", "m/s4",  0.5, Yc[:,_x, 4]),
             ("$y^{(4)}$", "m/s4",  0.5, Yc[:,_y, 4]),
             ("$z^{(4)}$", "m/s4",  0.5, Yc[:,_z, 4]),
             ("$\\psi^{(4)}$", "deg/s4",   0.5, np.rad2deg(Yc[:,_psi, 4])),
    ]
    figure, axes = ppu.plot_in_grid(time, plots, 4, figure, axes, window_title)
    return figure, axes

# /home/poine/work/two_d_guidance/two_d_guidance/path_factory.py
def plot2d(time, Yc, figure=None, window_title="Flat Output Trajectory"):
    figure = ppu.prepare_fig(figure, window_title, (20.48, 10.24))
    ax = plt.gca()
    points = Yc[:, _x:_z, 0].reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(0, len(points))
    lc = matplotlib.collections.LineCollection(segments, cmap='jet', norm=norm)
    lc.set_array(np.arange(len(points)))
    lc.set_linewidth(2)
    line = ax.add_collection(lc)
    figure.colorbar(line, ax=ax)
    ax.set_xlim(points[:,0, _x].min(), points[:,0, _x].max())
    ax.set_ylim(points[:,0, _y].min(), points[:,0, _y].max())
    ax.set_aspect('equal'); plt.title('2D')
    
