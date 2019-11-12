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

    def get_wind(self, pos, t):
        return [0, 0, 0]
