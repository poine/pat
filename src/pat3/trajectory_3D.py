import sys, os, math, numpy as np
import logging
import pdb

#import matplotlib.pyplot as plt

import pat3.utils as p3_u
#import pat3.plot_utils as p3_pu
#import pat3.algebra as p3_alg

class RefTraj:
    def __init__(self):
        pass

    def get_points(self):
        return []

class CircleRefTraj():
    def __init__(self, c=(0, 0, 0), r=40.):
        self.c, self.r = np.asarray(c), r
        self.alpha0, self.dalpha = 0, 2*np.pi

    def _zfun(self, alpha): return  np.zeros_like(alpha)
    #def _zdfun(self, alpha): return  np.zeros_like(alpha)
    def _zdfun(self, alpha, eps=1e-2, zdm=10.):
        return np.clip((self._zfun(alpha+eps)-self._zfun(alpha-eps))/(2*eps), -zdm, zdm)

    def offset(self, delta_center):
        self.c += delta_center
    def set_center(self, _center):
        self.c = np.asarray(_center)
    
    def get_points(self):
        alphas = np.linspace(self.alpha0, self.alpha0+self.dalpha, int(360*self.dalpha/np.pi/2))
        xs = self.c[0] + self.r*np.cos(alphas)
        ys = self.c[1] + self.r*np.sin(alphas)
        zs = self.c[2] + self._zfun(alphas)# - 5*np.ones_like(alphas)#
        return np.vstack((xs, ys, zs)).T
        
    def get_alpha(self, p): # FIXME alpha0
        cp = p-self.c
        alpha = np.arctan2(cp[1], cp[0])
        #print('{} {:.1f}'.format(alpha, np.deg2rad(alpha)))#looks good
        return alpha

    def _pt(self, alpha):
        return self.c+[self.r*np.cos(alpha), self.r*np.sin(alpha), self._zfun(alpha)]
    
    def get_point_ahead(self, p0, lookahead_dist):
        alpha0 = self.get_alpha(p0)
        alpha1 = alpha0 + lookahead_dist/self.r
        p1 = self._pt(alpha1) #self.c+[self.r*np.cos(alpha1), self.r*np.sin(alpha1), self._zfun(alpha1)]
        return p1

    def get_proj_and_carrot(self, p0, lookahead_dist, _v=17.):
        alpha0 = self.get_alpha(p0)
        alpha1 = alpha0 + lookahead_dist/self.r
        p1, p2 = self._pt(alpha0), self._pt(alpha1)
        zd = self._zdfun(alpha0)*_v/self.r  #
        return p1, p2, zd
        

class BankedCircleRefTraj(CircleRefTraj):
    def __init__(self, c=(0, 0, 0), r=40.):
        CircleRefTraj.__init__(self, c, r)

    def _zfun(self, alpha): return  20*np.sin(alpha+3*np.pi/2) 

class ZStepCircleRefTraj(CircleRefTraj):
    def __init__(self, c=(0, 0, 0), r=40.):
        CircleRefTraj.__init__(self, c, r)

    def _zfun(self, alpha):
        def __zfun(alpha, _a=2.):
            alpha = p3_u.norm_0_2pi(alpha)
            return _a if alpha < np.pi else -_a
        if type(alpha) is np.float64: return __zfun(alpha)
        ret = np.ones_like(alpha)
        for i, alpha in enumerate(alpha): ret[i] = __zfun(alpha)
        return  ret

class ZdStepCircleRefTraj(CircleRefTraj):
    def __init__(self, c=(0, 0, 0), r=40.):
        CircleRefTraj.__init__(self, c, r)

    def _zfun(self, alpha, _a=10.):
        def __zfun(alpha):
            alpha = p3_u.norm_0_2pi(alpha)
            if alpha < np.pi/2: _h = alpha*_a
            elif alpha < np.pi: _h = _a*np.pi/2
            elif alpha < 3*np.pi/2: _h = (3*np.pi/2-alpha)*_a
            else: _h=0 
            return _h#alpha*_a if alpha < np.pi else (2*np.pi-alpha)*_a
        if type(alpha) is np.float64: return __zfun(alpha)
        ret = np.ones_like(alpha)
        for i, alpha in enumerate(alpha):
            ret[i] = __zfun(alpha)
        return  ret
        
class LineRefTraj():
    def __init__(self, p1=(0, 0, 0), p2=(50, 0, 0)):
        self.p1, self.p2 = np.asarray(p1), np.asarray(p2)
        self.p1p2 = self.p2-self.p1
        self.dist =  np.linalg.norm(self.p1p2)
        self.n = self.p1p2/self.dist
        
    def get_points(self, res=0.1):
        n_pts = self.dist/res
        return np.stack([np.linspace(self.p1[i], self.p2[i], n_pts) for i in range(3)], axis=-1)
    
    def get_point_ahead(self, p, lookahead_dist):
        p1p = p-self.p1
        _np1p = np.linalg.norm(p1p)
        foo = max(0, np.dot(self.n, p1p)) # don't aim before starting point
        p0 = self.p1 + self.n*(foo + lookahead_dist)
        if np.dot(self.n, p1p) > self.dist: p0 = self.p2 + self.n*lookahead_dist
        return p0

    def finished(self, p, lookahead_dist):
        p1p = p-self.p1
        return  np.dot(self.n, p1p) > self.dist

    def dist_to_end(self, p):
        return self.dist - np.dot(self.n, p-self.p1)

        
    
class SpiralRefTraj:
    def __init__(self, c=(0, 0, 10), r=30.):
        pass


class CompositeRefTraj:
    def __init__(self, sub_trajs):
        self.sub_trajs = sub_trajs
        self.cur_traj_id = 0

    def get_points(self, res=0.1):
        return np.concatenate([s.get_points() for s in self.sub_trajs])

    def get_point_ahead(self, p, lookahead_dist):
        dte = self.sub_trajs[self.cur_traj_id].dist_to_end(p)
        if dte < lookahead_dist:
            next_traj_id = (self.cur_traj_id+1) % len(self.sub_trajs)
            print(self.cur_traj_id, next_traj_id)
            if dte <= 0:
                self.cur_traj_id = next_traj_id
                
            return self.sub_trajs[next_traj_id].get_point_ahead(p, lookahead_dist-dte)
        else:
            return self.sub_trajs[self.cur_traj_id].get_point_ahead(p, lookahead_dist)


class SquareRefTraj(CompositeRefTraj):
    def __init__(self):
        CompositeRefTraj.__init__(self, [LineRefTraj(p1=( 50,  50, 0), p2=( 50, -50, 0)),
                                         LineRefTraj(p1=( 50, -50, 0), p2=(-50, -50, 0)),
                                         LineRefTraj(p1=(-50, -50, 0), p2=(-50,  50, 0)),
                                         LineRefTraj(p1=(-50,  50, 0), p2=( 50,  50, 0))])

class OctogonRefTraj(CompositeRefTraj):
    def __init__(self, a=50):
        lines = []
        for i in range(8):
            a1, a2 = i*np.pi/3, (i+1)*np.pi/3
            lines.append(LineRefTraj(p1=( a*np.sin(a1), a*np.cos(a1), 0),
                                     p2=( a*np.sin(a2), a*np.cos(a2), 0)))
        CompositeRefTraj.__init__(self, lines)        


