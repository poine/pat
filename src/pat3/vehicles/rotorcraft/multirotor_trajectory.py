#-*- coding: utf-8 -*-

'''
  (Output) Trajectories for a quadrotor
  Y = [x, y, z, psi] plus their 4 time derivatives
  Y constitutes a flat output and can hence be used to compute the state and input
  of the quadrotor (see vehicle/rotorcraft/multirotor_control.py)
'''
import math, numpy as np, matplotlib.pyplot as plt, matplotlib
from abc import ABC, abstractmethod

import pat3.utils as p3_u, pat3.plot_utils as p3_pu, pat3.algebra as p3_al
import pat3.trajectory_1D as p3_t1D

_x, _y, _z, _psi, _ylen = range(5) # components
_nder = 5 # we consider 4 time derivatives, this should be a parameter...

class Trajectory(ABC):
    _x, _y, _z, _psi, _ylen = range(5)
    _nder = 5
    def is_space_indexed(self): return False  # space indexed trajectories return True
    def has_dyn_ctl_pts(self): return False   # space_indexed trajectories might have them  
    def get_dyn_ctl_pts(self): return None    #  
    def set_dyn_ctl_pts(self, pts): pass      #
    def has_waypoints(self): return False     # geometry control points (why not?)
    def get_waypoints(self): return None      # 
    def set_waypoints(selp, pts): pass        #
    def set_t0(self, t0): pass                # set start time (useful for sequencing)
    @abstractmethod
    def get(self, t): pass
    
class Cst(Trajectory):
    def __init__(self, Y00=None, duration=1.):
        self.duration = duration
        self.Y = np.zeros((_ylen, _nder))
        if Y00 is not None: self.Y[:,0] = Y00

    def get(self, t):
        return self.Y

class Circle(Trajectory):
    def __init__(self, c=[0, 0, 0], r=1., v=2., alpha0=0, dalpha=2*np.pi, zt=None, psit=None):
        self.c, self.r, self.v = np.asarray(c), r, v # center, radius, velocity
        self.alpha0, self.dalpha = alpha0, dalpha    # start angle, angle span
        self.omega = self.v/self.r                   # angular velocity
        self.length = dalpha*np.abs(r)               
        self.duration = self.length/v
        self.zt = zt if zt is not None else p3_t1D.CstOne(c[_z])
        self.psit = psit or p3_t1D.AffineOne(self.omega, self.alpha0+np.sign(self.r)*np.pi/2) # forward facing (enu?)
        #self.psit = psit or p3_t1D.AffineOne(self.omega, self.alpha0)
        self.t0 = 0.
        self.omp = [self.omega**_i for _i in range(_nder)] # precompute angular velocity powers
        
    def set_t0(self, t0): self.t0 = t0
     
    def get(self, t):
        dt = t-self.t0
        alpha =  self.alpha0 + self.omega*(dt)
        rca, rsa = np.abs(self.r)*np.cos(alpha), np.abs(self.r)*np.sin(alpha) 
        Yc = np.zeros((_nder,_ylen))
        ## xy
        Yc[0,:_z] = self.c[:_z] + [rca, rsa]
        Yc[1,:_z] = [-self.omp[1]*rsa,  self.omp[1]*rca]
        Yc[2,:_z] = [-self.omp[2]*rca, -self.omp[2]*rsa]
        Yc[3,:_z] = [ self.omp[3]*rsa, -self.omp[3]*rca]
        Yc[4,:_z] = [ self.omp[4]*rca,  self.omp[4]*rsa]
        ## z
        Yc[:,_z] = self.zt.get(dt)
        ## psi
        Yc[:,_psi] =  self.psit.get(dt)
        Yc[0,_psi]  = p3_al.norm_mpi_pi(Yc[0,_psi])
        return Yc.T

class Line(Trajectory):
    def __init__(self, p1, p2, v=2., psi=None):
        self.p1, self.p2, self.v = np.asarray(p1), np.asarray(p2), v # start, end and velocity 
        dep = self.p2-self.p1
        self.length = np.linalg.norm(dep)   # length
        self.un = dep/self.length           # unit vector
        self.psi = psi if psi is not None else np.arctan2(self.un[1], self.un[0])
        self.duration = self.length/self.v  # duration
        self.t0 = 0.

    def set_t0(self, t0): self.t0 = t0
        
    def get(self, t):
        Yc = np.zeros((_nder, _ylen))
        Yc[0,:3] = self.p1 + self.un*self.v*(t-self.t0)
        Yc[1,:3] =           self.un*self.v
        Yc[0,3] = self.psi
        return Yc.T


class SmoothLine:
    def __init__(self, Y00=[0, 0, 0, 0], Y10=[1, 0, 0, 0], duration=1.):
        self.duration = duration
        Y0 = np.zeros((_ylen, _nder))
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
        
    def set_t0(self, t0):
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
            s.set_t0(st)
        self.t0 = 0.   

    def set_t0(self, t0):
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
            self.refs[i].set_t0(self.Ysp[0,i] if Y0 is None else Y0[i])

        self.Ys = np.zeros_like(self.Ysp)
        for i in range(1, len(self.time)):
            for ycmp in range(4):
                if ycmp==3:  # wrap pesky heading angle
                    err = self.Ys[i-1, ycmp, 0] - self.Ysp[i, ycmp, 0]
                    if err > np.pi: self.Ysp[i, ycmp, 0] += 2*np.pi
                    elif err < -np.pi: self.Ysp[i, ycmp, 0] -= 2*np.pi
                self.Ys[i, ycmp] = self.refs[ycmp].run(dt, self.Ysp[i, ycmp, 0])

        
    def set_t0(self, t0): pass
        
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
    figure = p3_pu.prepare_fig(figure, window_title, (20.48, 10.24))
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
    figure, axes = p3_pu.plot_in_grid(time, plots, 4, figure, axes, window_title)
    return figure, axes

# /home/poine/work/two_d_guidance/two_d_guidance/path_factory.py
def plot2d(time, Yc, figure=None, window_title="Flat Output Trajectory"):
    figure = p3_pu.prepare_fig(figure, window_title, (20.48, 10.24))
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
    
