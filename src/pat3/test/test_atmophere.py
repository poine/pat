#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#

import numpy as np, matplotlib.pyplot as plt
import matplotlib, matplotlib.pyplot as plt, mpl_toolkits.mplot3d
from matplotlib import cm


import pat3.plot_utils as p3_plu, pat3.atmosphere as p3_atm


def plot_hslice(atm, h0=0.):
    
    p3_plu.plot_slice_wind_ne(atm, n0=-100, n1=100, dn=5., e0=-100., e1=100, de=5, h0=h0, t0=0.,
                              show_color_bar=True)
    plt.savefig('/tmp/atm2.png')

def plot_vslice(atm, h0=0, h1=250, n0=-200, n1=200, show_quiver=False):
    #atm = p3_atm.AtmosphereVgradient()
    #atm = p3_atm.AtmosphereThermal()
    #atm = p3_atm.AtmosphereThermal1()
    #atm =  p3_atm.AtmosphereWharington()
    #p3_plu.plot_slice_wind(atm, xmax=100, dx=2.5, h1=200)
    #atm = p3_atm.AtmosphereRidge()
    #p3_plu.plot_slice_wind (atm, xmax=40, dx=0.5, h0=-40, h1=40, dh=0.5)
    #p3_plu.plot_slice_wind2(atm, xmax=40, dx=1., h0=-40, h1=40, dh=1.)
    #p3_plu.plot_3D_wind(atm) # does not work
    p3_plu.plot_slice_wind_nu(atm, n0=n0, n1=n1, dn=5., e0=0., h0=h0, h1=h1, dh=2., zdir=-1.,
                              show_quiver=show_quiver, show_color_bar=True,
                              figure=None, ax=None)


def display_all_models(h0=0.):
    #atm1 = p3_atm.AtmosphereVgradient()
    atm1 =  p3_atm.AtmosphereWharington()
    atm2 =  p3_atm.AtmosphereGedeon()
    atm3 =  p3_atm.AtmosphereAllen()
    #atm3 = p3_atm.AtmosphereRidge()
    #atm4 = p3_atm.AtmosphereThermal()
    atms, names = [atm1, atm2, atm3], ['Wharington', 'Gedeon', 'Allen']
    for atm, name in zip(atms, names):
        #p3_plu.plot_slice_wind_ne(atm, n0=-100, n1=100, dn=5., e0=-100., e1=100, de=5, h0=h0, t0=0.,
        #                          show_color_bar=True) 
        p3_plu.plot_slice_wind_nu(atm, n0=-100, n1=100, dn=5., e0=0., h0=-10, h1=180, dh=2., zdir=-1.,
                                  show_quiver=True, show_color_bar=True,
                                  title=name,
                                  figure=None, ax=None)

def test_ridge():
    atm = p3_atm.AtmosphereRidge()
    p3_plu.plot_slice_wind_nu(atm, n0=-60, n1=100, dn=10., e0=0., h0=0, h1=100, dh=10., zdir=-1.,
                              show_quiver=True, show_color_bar=True,
                              figure=None, ax=None)
    p3_plu.plot_slice_wind_eu(atm, n0=0, e0=-50, e1=50, de=10., h0=0, h1=100, dh=10., zdir=-1.,
                              show_quiver=True, show_color_bar=True,
                              figure=None, ax=None)

def test_wharington():
    atm = p3_atm.AtmosphereWharington()
    p3_plu.plot_slice_wind_nu(atm, n0=-60, n1=100, dn=10., e0=0., h0=0, h1=100, dh=10., zdir=-1.,
                              show_quiver=True, show_color_bar=True,
                              figure=None, ax=None)
    # p3_plu.plot_slice_wind_eu(atm, n0=0, e0=-50, e1=50, de=10., h0=0, h1=100, dh=10., zdir=-1.,
    #                           show_quiver=True, show_color_bar=True,
    #                           figure=None, ax=None)


# FIXME: this shows that grids need to be interpolated
def test_grid():
    nc_f = '/home/poine/work/glider_experiments/data/extr_IHODC.1.RK4DI.007.nc'
    atm = p3_atm.AtmosphereNC(nc_f)
    print(f'range n:{atm.x0}   {atm.x1} resolution: {atm.dx/(len(atm.ni)-1)}')
    print(f'range e:{atm.y0}   {atm.y1} resolution: {atm.dy/(len(atm.nj)-1)}')
    print(f'range d:{atm.z0}   {atm.z1} resolution: {atm.dz/(len(atm.level)-1)}')

    # North oriented line
    ns, e, h = np.arange(0, 200), 25, 100
    wzs =  [atm.get_wind_ned1([n, e, -h], 0)[2] for n in ns]
    wzs2 =  [atm.get_wind_ned3([n, e, -h], 0)[2] for n in ns]
    plt.plot(ns, wzs, '.', label='1')
    plt.plot(ns, wzs2, '.', label='interp')
    p3_plu.decorate(plt.gca(), f'east: {e} height:{h}', 'north in m', 'vz in m/s', legend=True)
    if 0:
        p3_plu.plot_slice_wind_nu(atm, ns[0], ns[-1], dn=5., e0=e, h0=0, h1=100, dh=2., zdir=-1.,
                                  show_quiver=True, show_color_bar=True,
                                  figure=None, ax=None)
    plt.figure()
    # East oriented line
    n, es, h = 25, np.arange(0, 200), 100
    wzs =  [atm.get_wind_ned1([n, e, -h], 0)[2] for e in es]
    wzs2 =  [atm.get_wind_ned3([n, e, -h], 0)[2] for e in es]
    plt.plot(es, wzs, '.', label='1')
    plt.plot(es, wzs2, '.', label='interp')
    p3_plu.decorate(plt.gca(), f'north: {n} height:{h}', 'east in m', 'vz in m/s', legend=True)

    plt.figure()
    # up oriented line
    n, e, hs = 25, 25, np.arange(0, 300)
    wzs =  [atm.get_wind_ned1([n, e, -h], 0)[2] for h in hs]
    wzs2 =  [atm.get_wind_ned3([n, e, -h], 0)[2] for h in hs]
    plt.plot(hs, wzs, '.', label='1')
    plt.plot(hs, wzs2, '.', label='interp')
    p3_plu.decorate(plt.gca(), f'north: {n} east:{e}', 'h in m', 'vz in m/s', legend=True)



        
def main():
    #atm =  p3_atm.AtmosphereWharington()
    #atm =  p3_atm.AtmosphereWharingtonArray()
    #atm =  p3_atm.AtmosphereAllen()
    #plot_vslice(atm, h0=-2000, h1=250)
    #atm = p3_atm.AtmosphereRidge()
    #plot_vslice(atm, h0=-10, h1=100, n0=-50, n1=150, show_quiver=True)
    #plot_hslice(atm)
    #display_all_models()
    #test_ridge()
    #test_wharington()
    test_grid()
    plt.show()



    
if __name__ == "__main__":
    main()
  
    

