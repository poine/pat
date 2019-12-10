#!/usr/bin/env python
#-*- coding: utf-8 -*-
#
# Copyright 2013-2015 Antoine Drouin (poinix@gmail.com)
#
# This file is part of PAT.
#
#    PAT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PAT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PAT.  If not, see <http://www.gnu.org/licenses/>.
#

import pdb

"""
  atmophère, atmosphère, est-ce-que j'ai une gueule d'atmosphère?
"""

import math, numpy as np, matplotlib.pyplot as plt

def mach_of_va(va, T, k=1.4, Rs=287.05): return va/math.sqrt(k*Rs*T)

def get_rho(h):
    """
    Get air density at the given altitude (in m)
    Warning, the two parts of that function do not
    join not well!
    """
    if h<= 11000:
        return 1.225 * math.pow((1-6.5e-3*h/288.15), 4.2557)
    else:
        return 0.36 * math.exp(-1.17e-4*(h-11000))

'''
International Standard Atmosphere Model
see: http://en.wikipedia.org/wiki/International_Standard_Atmosphere
'''
_name, _h0, _z0, _a, _T0, _p0 = np.arange(6)
#      name         h(m)    z(km)    a(K/m)    T0(K)    p0(Pa)
param = \
[['Troposphere',       0,   0.0,    -6.5e-3,   288.15,   101325],
 ['Tropopause',    11000,  11.019,   0.,       216.65,    22632],
 ['Stratosphere',  20000,  20.063,   1.0e-3,   216.65,     5474.9],
 ['Stratosphere',  32000,  32.162,   2.8e-3,   228.65,      868.02],
 ['Stratopause',   47000,  47.350,   0.0e-3,   270.65,      110.91],
 ['Mesosphere',    51000,  51.413,  -2.8e-3,   270.65,    66.939],
 ['Mesosphere',    71000,  71.802,  -2.0e-3,   214.65,     3.9564],
 ['Mesopause',     84852,  86.000,   0.,       186.87,     0.3734]]

def isa(h):
    layer = 0
    while param[layer][_h0] < h: layer+=1
    name, h0, z0, a, T0, p0 = param[layer-1]
    dh = h - h0
    T = T0 + a*dh
    g=9.81; R=287.0
    if a != 0.:
        p = p0*math.pow(T/T0, -g/a/R)
    else:
        p = p0*math.exp(-g/R/T0*dh)
    rho = p/R/T
    return p, rho, T

def decorate(ax, title=None, xlab=None, ylab=None, legend=None):
    ax.xaxis.grid(color='k', linestyle='-', linewidth=0.2)
    ax.yaxis.grid(color='k', linestyle='-', linewidth=0.2)
    if xlab:
        ax.xaxis.set_label_text(xlab)
    if ylab:
        ax.yaxis.set_label_text(ylab)
    if title:
        ax.set_title(title, {'color'    : 'k', 'fontsize'   : 20 })
    if legend != None:
        ax.legend(legend, loc='best')


def plot(h0=1, h1=84000):
    h = np.linspace(h0, h1, 1000)
    v = np.array(map(isa, h))
    v1 = np.array(map(get_rho, h))
    print(v1)
    #pdb.set_trace()
    ax = plt.subplot(1, 3, 1)
    plt.plot(v[:,2], h)
    decorate(ax, 'Temperature', 'K', 'm')
    ax = plt.subplot(1, 3, 2)
    plt.plot(v[:,1], h)
    plt.plot(v1, h, 'r')
    decorate(ax, 'Density', 'Kg/m3', 'm')
    ax = plt.subplot(1, 3, 3)
    plt.plot(v[:,0], h)
    decorate(ax, 'Pressure', 'Pa', 'm')
    plt.show()

if __name__ == "__main__":
    h = 1
    print("isa at {}m : {} Pa {} kg/m3 {} K".format(h, *isa(h)))
    #plot(h0=1, h1=84000)
    plot(h0=1, h1=11000)




class Atmosphere:

    def sample(self, x_r=np.arange(-60., 60., 5), y_r=np.arange(-60., 60., 5), z_r=np.arange(0., -150, -10), t=0.):
        x, y, z = np.meshgrid(x_r, y_r, z_r)
        wx, wy, wz = np.meshgrid(x_r, y_r, z_r)
        nx, ny, nz = x.shape
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    pos = [x[ix, iy, iz], y[ix, iy, iz], z[ix, iy, iz]]
                    wx[ix, iy, iz], wy[ix, iy, iz], wz[ix, iy, iz] = self.get_wind(pos, t)
        return x, y, z, wx, wy, wz

    
    def get_wind(self, pos_ned, t):
        return [0, 0, 0]
        #return [0., 0, 0.5] if pos[0] > 0 else [0, 0, -0.5]
        #d = np.linalg.norm(pos-[10, 10, 0])
        #return [0, 0, np.sin(0.1*d)]
        #return [0., 0, 0.5] if pos[0] > 0 else [0, 0, -0.5]
    


class AtmosphereCalm(Atmosphere):
    def get_wind(self, pos_ned, t):
        return [0, 0, 0] 

class AtmosphereCstWind(Atmosphere):
    def __init__(self, v=[0, 0, 0]):
        self.v = np.asarray(v)
    def get_wind(self, pos, t):
        #p = 120.
        #return self.v if math.fmod(t, p) > p/3 else [0, 0, 0]
        return self.v

class AtmosphereVgradient(Atmosphere):
    def get_wind(self, pos, t):
        wmin, wmax, hmax = 0., 5., 10. # 
        return [max(wmin, min(wmax/hmax*pos[2], wmax)), 0, 0]


def thermal_model_gedeon(x, y, z, R_t=10, W_max=1):
    r = np.sqrt(x**2+y**2)
    top = (r/R_t)**2
    w = W_max*np.exp(-top)*(1-top)
    return w

class AtmosphereThermal(Atmosphere):
    def get_wind(self, pos_ned, t, wstar=256, wgain=1., rgain=1, zi=1000):
        x,y,z = pos_ned
        z=-z
        w = thermal_model_gedeon(x, y, z, R_t=30, W_max=2)
        return [0, 0, -w]


class AtmosphereThermal1(Atmosphere):
    def __init__(self):
        self.c = np.array([0, 0, 0])
        self.zi = 2000
        self.wstar = 256.

    def set_params(self, xc, yc, zi, wstar):
        self.c = np.asarray([xc, yc, 0])
        self.zi, self.wstar = zi, wstar

    def get_params(self): return self.c[0], self.c[1], self.zi, self.wstar  
        
    def get_wind(self, pos_ned, t, wgain=1., rgain=1):
        x,y,z = pos_ned - self.c
        z=z+150
        #CALCULATE AVERAGE UPDRAFT SIZE
        zzi=z/self.zi
        rbar=(.102*zzi**(1./3))*(1-(.25*zzi))*self.zi
        #CALCULATE AVERAGE UPDRAFT STRENGTH
        wtbar=np.power(zzi,1./3)*(1-1.1*zzi)*self.wstar
        
        #CALCULATE INNER AND OUTER RADIUS OF ROTATED TRAPEZOID UPDRAFT
        r2=rbar*rgain #multiply by random perturbation gain
        if r2<10: r2=10
        if r2<600:
            r1r2=.0011*r2+.14
        else:
            r1r2=.8
        r1=r1r2*r2
        #limit small updrafts to 20m diameter
        #MULTIPLY AVERAGE UPDRAFT STRENGTH BY WGAIN FOR THIS UPDRAFT
        wt=wtbar*wgain #add random perturbation
        #CALCULATE STRENGTH AT CENTER OF ROTATED TRAPEZOID UPDRAFT
        wc=(3*wt*((r2**3)-(r2**2)*r1)) / ((r2**3)-(r1**3))
        #print wc
        wc/=300.
        R_ave = (r1+r2)/2.
        w = thermal_model_gedeon(x, y, z, R_t=R_ave, W_max=wc)
        return [0, 0, -w]

class AtmosphereThermalMoving(AtmosphereThermal1):

    def get_wind(self, pos_ned, t):
        self.c = np.array([0, 50*np.sin(0.2*t), 0])
        return AtmosphereThermal1.get_wind(self, pos_ned, t)

class AtmosphereThermalMulti(Atmosphere):
    def __init__(self):
        self.thermals = [AtmosphereThermal1() for i in range(2)]
        zi, wstar =850., 256.
        self.thermals[0].set_params(0,  55, zi, wstar)
        self.thermals[1].set_params(0, -55, zi, wstar)
        
    def set_params(self, xc, yc, zi, wstar, idx=0):
        print('set params:', xc, yc, zi, wstar, idx)
        self.thermals[idx].set_params(xc, yc, zi, wstar)

    def get_params(self, idx=0):
        return self.thermals[idx].get_params()
    
    def get_wind(self, pos_ned, t): 
        #self.thermals[0].c = np.array([0, 55+10*np.sin(0.02*t), 0])
        winds = [_t.get_wind(pos_ned, t) for _t in self.thermals]
        return np.sum(winds, axis=0)



class AtmosphereRidge(Atmosphere):
    # wind over a cylinder obstacle.
    # see: Langellan, long distance/duration trajectory optimization for small uavs
    def __init__(self):
        self.R = 30            # cylinder radius
        self.winf = 2.         # 
        self.R2 = self.R**2
        self.c = np.array([40, 0, 15])
        
    def set_params(self, *args): pass


    def get_wind(self, pos_ned, t): 
        dpos = pos_ned - self.c
        dpos[1]=0
        r = np.linalg.norm(dpos)
        if r < self.R: return [0, 0, 0]
        eta = -np.arctan2(dpos[2], dpos[0])
        ceta, seta = np.cos(eta), np.sin(eta)
        R2ovr2 = self.R2/r**2
        wx = self.winf*(1-R2ovr2*(ceta**2-seta**2))
        wz = 2*self.winf*R2ovr2*ceta*seta
        return [wx, 0, wz]
