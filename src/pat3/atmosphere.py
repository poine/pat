#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#
# Copyright 2013-2021 Antoine Drouin (poinix@gmail.com)
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
  atmosphere related stuff
   -isa model
   -misc wind distributions
"""

import math, numpy as np, matplotlib.pyplot as plt
import scipy.interpolate
try:
    import netCDF4
except ImportError:
    print('pat3.atmosphere: netcfd4 not available on your system. You won\'t be able to load grids')

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

    ## FIXME: remove this function, call get_wind_ned instead
    def get_wind(self, pos_ned, t): return self.get_wind_ned(pos_ned, t)

    def get_wind_ned(self, pos_ned, t): return [0, 0, 0]
    def get_wind_enu(self, pos_enu, t):
        pos_ned = [pos_enu[1], pos_enu[0], -pos_enu[2]]
        wind_ned = self.get_wind_ned(self, pos_ned, t)
        wind_enu = [wind_ned[1], wind_ned[0], -wind_ned[2]]
        return wind_enu


class AtmosphereCalm(Atmosphere):
    def get_wind_ned(self, pos_ned, t):
        return [0, 0, 0] 

class AtmosphereCstWind(Atmosphere):
    def __init__(self, v=[0, 0, 0]):
        self.v = np.asarray(v)
    def get_wind_ned(self, pos, t):
        return self.v

class AtmosphereSinetWind(AtmosphereCstWind):
    def get_wind_ned(self, pos, t, omega=1.):
        return np.sin(omega*t)*self.v

class AtmosphereSinedWind(AtmosphereSinetWind):
    def get_wind_ned(self, pos, t, omega=0.1):
        return np.sin(omega*np.linalg.norm(pos))*self.v

class AtmosphereVgradient(Atmosphere):
    def __init__(self, w0=-2, w1=5, h0=20, h1=60):
        self.w0, self.w1, self.h0, self.h1 = w0, w1, h0, h1
        self.dw, self.dh = w1-w0, h1-h0
        
    def get_wind_ned(self, pos, t):
        dh = -pos[2]-self.h0
        wx = np.clip(self.w0+dh/self.dh*self.dw, self.w0, self.w1)
        return [-wx, 0, 0]


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

# class AtmosphereThermalMulti(Atmosphere):
#     def __init__(self):
#         self.thermals = [AtmosphereThermal1() for i in range(2)]
#         zi, wstar =850., 256.
#         self.thermals[0].set_params(0,  55, zi, wstar)
#         self.thermals[1].set_params(0, -55, zi, wstar)
        
#     def set_params(self, xc, yc, zi, wstar, idx=0):
#         print('set params:', xc, yc, zi, wstar, idx)
#         self.thermals[idx].set_params(xc, yc, zi, wstar)

#     def get_params(self, idx=0):
#         return self.thermals[idx].get_params()
    
#     def get_wind(self, pos_ned, t): 
#         #self.thermals[0].c = np.array([0, 55+10*np.sin(0.02*t), 0])
#         winds = [_t.get_wind(pos_ned, t) for _t in self.thermals]
#         return np.sum(winds, axis=0)



#
# Shear Model
#
class AtmosphereShearX(Atmosphere):
    # Wind shear, only in X.
    # Adapted from Drela's DSOpt
    def __init__(self, wind1=15.0, wind2=-2.0, xlayer=30.0, zlayer=40.0,txcon=0.11):
        self.wind1  = wind1    #   wind speed above shear layer  [m/s] ,   V[m/s] = V[mph] * 0.44694
        self.wind2  = wind2    #   wind speed below shear layer  [m/s]
        self.xlayer = xlayer   #   shear layer origin distance from orbit center     [m]
        self.zlayer = zlayer   #   shear layer centerline height above orbit center  [m]
        self.txcon  = txcon    # 0.11  shear layer spreading constant (0.11 is from mixing-layer theory)
        # shear layer spreading rate  d(thickness)/dx
        self.tlayerx = 2.0*self.txcon * (self.wind1-self.wind2)/(self.wind1+self.wind2) # wind1 should not be = to -wind2 !!!
        
    def set_params(self, *args): pass


    def get_wind_ned(self, pos_ned, t):
        ''' Outputs wind vector : wind = [windx,windy,windz]'''
        xdel = pos_ned[0] - self.xlayer
        zdel = -pos_ned[2] - self.zlayer

        tlayer = -xdel*self.tlayerx
        tlayer = max( tlayer , 1.0e-5 )

        znorm = zdel / tlayer
        znorm = max( -0.5 , min( 0.5 , znorm ) )

        frac = 0.5 + 0.5*np.sin(np.pi*znorm)
        wx = -self.wind1*frac - self.wind2*(1.0-frac) 
        wy = 0.
        wz = 0.
        return [wx,wy,wz]

#
# Gaussian model
#
class AtmosphereWharington(Atmosphere):
    
    def __init__(self, center=None, radius=50, strength=-2, cst=[0, 0, 0]):
        self.center = np.asarray(center) if center is not None else np.array([0, 0, 0])
        self.radius = radius
        self.strength = strength
        self.cst = cst
        self.r2 = self.radius**2
        #self.x0, self.x1, self.y0, self.y1, self.z0, self.z1 = -100, 100, -100, 100, 0, 100
        
    def set_params(self, *args): pass

    def get_wind_ned(self, pos_ned, t): 
        dpos = pos_ned - self.center
        r2 = dpos[0]**2+dpos[1]**2
        wz = self.strength*np.exp(-r2/self.r2)
        return np.array([0, 0, wz])+self.cst

class AtmosphereWharingtonArray(Atmosphere):
    def __init__(self, centers, radiuses, strengths):
        self.thermals = [AtmosphereWharington(_c, _r, _s) for _c, _r, _s in zip(centers, radiuses, strengths)]

    def get_wind_ned(self, pos_ned, t): 
        winds = [_t.get_wind_ned(pos_ned, t) for _t in self.thermals]
        return np.sum(winds, axis=0)

    def set_params(self, *args): pass


class AtmosphereArray(Atmosphere):
    def __init__(self, atms):
        self.atms = atms
    
    def get_wind_ned(self, pos_ned, t): 
        winds = [_t.get_wind_ned(pos_ned, t) for _t in self.atms]
        return np.sum(winds, axis=0)

    def set_params(self, *args): pass
    
        
# class AtmosphereThermalMulti(AtmosphereWharingtonArray):
#     def __init__(self):
#         Cs, Rs, Ss = [[75, -20, 0], [-75, 20, 0]], [75, 50], [-0.5, -1.2]
#         #Cs, Rs, Ss = [[0, 0, 0]], [50], [-2]
#         AtmosphereWharingtonArray.__init__(self, Cs, Rs, Ss)

#
# Improved Gaussian model (with downdraft)
#
class AtmosphereGedeon(AtmosphereWharington):

    def get_wind_ned(self, pos_ned, t): 
        dpos = pos_ned - self.center
        r2 = dpos[0]**2+dpos[1]**2
        top = r2/self.r2
        wz =  self.strength*np.exp(-top)*(1-top)
        return np.array([0, 0, wz])

class AtmosphereGedeonArray(Atmosphere):
    def __init__(self, centers, radiuses, strengths):
        self.thermals = [AtmosphereGedeon(_c, _r, _s) for _c, _r, _s in zip(centers, radiuses, strengths)]

    def get_wind_ned(self, pos_ned, t): 
        winds = [_t.get_wind_ned(pos_ned, t) for _t in self.thermals]
        return np.sum(winds, axis=0)

    def set_params(self, *args): pass

class AtmosphereThermalMulti(AtmosphereGedeonArray):
     def __init__(self):
         Cs, Rs, Ss = [[75, -20, 0], [-75, 20, 0]], [75, 50], [-0.5, -1.2]
         AtmosphereGedeonArray.__init__(self, Cs, Rs, Ss)
    

# Wind over a cylinder obstacle.
# see: Langellan, long distance/duration trajectory optimization for small uavs
class AtmosphereRidge(Atmosphere):
    def __init__(self, winf=2.):
        self.R = 50            # cylinder radius (was 30)
        self.winf = winf       # 
        self.R2 = self.R**2
        self.c = np.array([40, 0, 15])
        
    def set_params(self, *args): pass

    def get_wind_ned(self, pos_ned, t): 
        dpos = pos_ned - self.c
        dpos[1]=0  # cylinder axis is y (east)
        r = np.linalg.norm(dpos)
        if r < self.R: return [0.,0.,0.]
        eta = -np.arctan2(dpos[2], dpos[0])
        ceta, seta = np.cos(eta), np.sin(eta)
        R2ovr2 = self.R2/r**2
        wx = self.winf*(1-R2ovr2*(ceta**2-seta**2))
        wz = 2*self.winf*R2ovr2*ceta*seta
        return [wx, 0, wz]


    
#
# Allen model
#
# This needs to be fixed, not sure i understand the equations :(
class AtmosphereAllen(Atmosphere):
    def __init__(self, center=None):
        self.center = np.asarray(center) if center is not None else np.array([0, 0, 0])
        self.zi = 2000
        self.wstar = 256.
        self.dz = 200.
        self.x0, self.x1, self.y0, self.y1, self.z0, self.z1 = -100, 100, -100, 100, 0, 100
       
    def get_wind(self, pos_ned, t, rgain=1., wgain=1.):
        x,y,z = pos_ned - self.center
        if z < -self.dz: return [0, 0, 0]
        z=z+self.dz
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

    

#from Scientific.Functions.Interpolation import InterpolatingFunction
from scipy import interpolate

# NetCDF4 grid
class AtmosphereNC(Atmosphere):
    def __init__(self, filename, center=None):
        self.center = np.asarray(center) if center is not None else np.array([0, 0, 0])
        self.nc_f = netCDF4.Dataset(filename, 'r')
        self.UT, self.VT, self.WT = [self.nc_f.variables[what][:] for what in ['UT', 'VT', 'WT']]
        # nc_attrs, nc_dims, nc_vars = analyse_nc(self.nc_f)
        self.ni, self.nj, self.level = [self.nc_f.variables[what][:] for what in ['ni', 'nj', 'level']]
        self.x0, self.x1 = self.ni[0], self.ni[-1]; self.dx = self.x1-self.x0
        self.y0, self.y1 = self.nj[0], self.nj[-1]; self.dy = self.y1-self.y0
        _id, _z = np.arange(len(self.level)), self.level
        self.level_of_z = scipy.interpolate.interp1d(_z, _id) # really?
        self.z0, self.z1 = self.level[0], self.level[-1]; self.dz=self.z1-self.z0

        self.interp_UT = scipy.interpolate.RegularGridInterpolator((self.level, self.nj, self.ni), self.UT.squeeze(), bounds_error=False, fill_value=0.)
        self.interp_VT = scipy.interpolate.RegularGridInterpolator((self.level, self.nj, self.ni), self.VT.squeeze(), bounds_error=False, fill_value=0.)
        self.interp_WT = scipy.interpolate.RegularGridInterpolator((self.level, self.nj, self.ni), self.WT.squeeze(), bounds_error=False, fill_value=0.)

        self.interp_modes = self.imod_none, self.imod_full_vz, self.imod_full_full = np.arange(3)
        self.interp_mode = self.imod_full_full
 
    
    def get_wind_ned(self, pos_ned, t):
        if self.interp_mode == self.imod_none:
            return self.get_wind_ned1(pos_ned, t)
        return self.get_wind_ned3(pos_ned, t)

    def get_wind_ned1(self, pos_ned, t):
        ni, nj, level = 0, 0, 0
        ni = int((pos_ned[0]-self.x0)/self.dx*len(self.ni))
        nj = int((pos_ned[1]-self.y0)/self.dy*len(self.nj))
        #ni = int((pos_ned[0]-self.x0)/self.dx*len(self.ni)+0.5)
        #nj = int((pos_ned[1]-self.y0)/self.dy*len(self.nj)+0.5)
        #pdb.set_trace()
        z = np.clip(-pos_ned[2], self.z0, self.z1)
        level = int(self.level_of_z(z))
        #print(z, level)
        #pdb.set_trace()
        _t=0 # TODO: fixme
        return [0, 0, self.WT[_t, level, nj, ni]]

    def get_wind_ned2(self, pos_ned, t):
        ni, nj, level = 0, 0, 0
        ni = (pos_ned[0]-self.x0)/self.dx*len(self.ni)
        nj = (pos_ned[1]-self.y0)/self.dy*len(self.nj)
        z = np.clip(-pos_ned[2], self.z0, self.z1)
        level = self.level_of_z(z)
        #print(z, level)
        _t=0 # TODO: add time...
        # no interpolation
        if 0:
            WT = self.WT[_t, int(level), int(nj), int(ni)]
        # interpolate north only
        if 0:
            ni1 = int(ni); ni2=ni1+1
            W1, W2 = self.WT[_t, int(level), int(nj), ni1], self.WT[_t, int(level), int(nj), ni2] 
            f = interpolate.interp1d([ni1, ni2], [W1, W2], assume_sorted=True)
            WT = float(f(ni))
        # interpolate north and east
        if 1:
            _xs, _ys = [int(ni), int(ni)+1], [int(nj), int(nj)+1]
            _zs = np.zeros((len(_xs), len(_ys)))
            for _ix, _x in enumerate(_xs):
                for _iy, _y in enumerate(_ys):
                    _zs[_ix, _iy] = self.WT[_t, int(level), _y, _x]
            _int = scipy.interpolate.RectBivariateSpline(_xs, _ys, _zs, kx=1, ky=1)
            WT = float(_int(ni, nj))
        #pdb.set_trace()
        return [0, 0, WT]

    def get_wind_ned3(self, pos_ned, t):
        hen = [-pos_ned[2], pos_ned[1], pos_ned[0]]
        if self.interp_mode == self.imod_full_full:
            UT = float(self.interp_UT(hen))
            VT = float(self.interp_VT(hen))
        else:
            UT, VT = 0, 0
        WT = float(self.interp_WT(hen))
        return [VT, UT, WT]

