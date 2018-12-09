#-*- coding: utf-8 -*-

import math, numpy as np, matplotlib.pyplot as plt, matplotlib
import pdb


'''
  Trajectories
  Y = [x, y, z, psi] plus their 4 time derivatives
'''
import pat3.plot_utils as ppu

_x, _y, _z, _psi = range(4)

class CstOne:
    def __init__(self, c=-1.):
        self.c = c
    def get(self, t):
        return np.array([self.c, 0, 0, 0, 0])

class AffineOne:
    def __init__(self, c1=-1., c2=0):
        self.c1, self.c2 = c1, c2
    def get(self, t):
        return np.array([self.c1*t+self.c2, self.c1, 0, 0, 0])
    
class SinOne:
    def __init__(self, c=0., a=1., om=1.):
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


    
class Circle:

    def __init__(self, c=[0, 0, 0], r=1., v=2., alpha0=0, dalpha=2*np.pi, zt=None, psit=None):
        self.c, self.r, self.v = np.asarray(c), r, v # center, radius, velocity
        self.alpha0, self.dalpha = alpha0, dalpha    # start angle, angle span
        self.omega = self.v/self.r
        self.length = dalpha*np.abs(r)
        self.duration = self.length/v
        self.zt = zt if zt is not None else CstOne(c[_z])
        self.psit = psit if psit is not None else AffineOne(self.omega, self.alpha0+np.sign(self.r)*np.pi/2)
        #self.psit = psit if psit is not None else AffineOne(self.omega, self.alpha0)
        self.t0 = 0.
        
    def reset(self, t0): self.t0 = t0
     
    def get(self, t):
        dt = t-self.t0
        alpha = self.omega*(dt) + self.alpha0
        rca, rsa = np.abs(self.r)*np.cos(alpha), np.abs(self.r)*np.sin(alpha) 
        Yc = np.zeros((5,4))
        ## xy
        Yc[0,:3] = self.c + [rca, rsa, 0]
        Yc[1,:3] = [-self.omega   *rsa,  self.omega   *rca, 0]
        Yc[2,:3] = [-self.omega**2*rca, -self.omega**2*rsa, 0]
        Yc[3,:3] = [ self.omega**3*rsa, -self.omega**3*rca, 0]
        Yc[4,:3] = [ self.omega**4*rca,  self.omega**4*rsa, 0]
        ## z
        Yc[:,_z] = self.zt.get(dt)
        ## psi
        Yc[:,_psi] = self.psit.get(dt)
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
        Yc = np.zeros((5,4))
        Yc[0,:3] = self.p1 + self.un*self.v*(t-self.t0)
        Yc[1,:3] =           self.un*self.v
        Yc[0,3] = self.psi
        return Yc.T
        
class Ellipse:
    pass

    
    


    
class CompositeTraj:
    def __init__(self, steps):
        self.steps = steps
        self.steps_dur = [s.duration for s in self.steps]
        self.steps_end = np.cumsum(self.steps_dur)
        self.duration = np.sum(self.steps_dur)
        for s, st in zip(self.steps[1:], self.steps_end):
            s.reset(st)
        self.t0 = 0.   

    def reset(self, t0): self.t0 = t0
        
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
        
    



def plot(time, Yc, figure=None, window_title="Flat Output Trajectory"):
    figure = ppu.prepare_fig(figure, window_title, (20.48, 10.24))
    #pdb.set_trace()
    plots = [("$x$",       "m",     0.5, Yc[:,_x, 0]),
             ("$y$",       "m",     0.5, Yc[:,_y, 0]),
             ("$z$",       "m",     0.5, Yc[:,_z, 0]),
             ("$\psi$",    "deg",   0.5, np.rad2deg(Yc[:,_psi, 0])),
             ("$x^{(1)}$", "m/s",   0.5, Yc[:,_x, 1]),
             ("$y^{(1)}$", "m/s",   0.5, Yc[:,_y, 1]),
             ("$z^{(1)}$", "m/s",   0.5, Yc[:,_z, 1]),
             ("$\psi^{(1)}$", "deg/s",   0.5, np.rad2deg(Yc[:,_psi, 1])),
             ("$x^{(2)}$", "m/s2",  0.5, Yc[:,_x, 2]),
             ("$y^{(2)}$", "m/s2",  0.5, Yc[:,_y, 2]),
             ("$z^{(2)}$", "m/s2",  0.5, Yc[:,_z, 2]),
             ("$\psi^{(2)}$", "deg/s2",   0.5, np.rad2deg(Yc[:,_psi, 2])),
             ("$x^{(3)}$", "m/s3",  0.5, Yc[:,_x, 3]),
             ("$y^{(3)}$", "m/s3",  0.5, Yc[:,_y, 3]),
             ("$z^{(3)}$", "m/s3",  0.5, Yc[:,_z, 3]),
             ("$\psi^{(3)}$", "deg/s3",   0.5, np.rad2deg(Yc[:,_psi, 3])),
             ("$x^{(4)}$", "m/s4",  0.5, Yc[:,_x, 4]),
             ("$y^{(4)}$", "m/s4",  0.5, Yc[:,_y, 4]),
             ("$z^{(4)}$", "m/s4",  0.5, Yc[:,_z, 4]),
             ("$\psi^{(4)}$", "deg/s4",   0.5, np.rad2deg(Yc[:,_psi, 4])),
    ]
    figure = ppu.plot_in_grid(time, plots, 4, figure, window_title)
    return figure

# /home/poine/work/two_d_guidance/two_d_guidance/path_factory.py
def plot3d(time, Yc, figure=None, window_title="Flat Output Trajectory"):
    figure = ppu.prepare_fig(figure, window_title, (20.48, 10.24))
    ax = plt.gca()
    points = Yc[:, _x:_z, 0].reshape(-1, 1, 2)
    #pdb.set_trace()
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
    
