#!/usr/bin/env python
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

def main():
    #atm =  p3_atm.AtmosphereWharington()
    #atm =  p3_atm.AtmosphereWharingtonArray()
    #atm =  p3_atm.AtmosphereAllen()
    #plot_vslice(atm, h0=-2000, h1=250)
    atm = p3_atm.AtmosphereRidge()
    plot_vslice(atm, h0=-10, h1=100, n0=-50, n1=150, show_quiver=True)
    #plot_hslice(atm)
    #display_all_models()
    plt.show()



    
if __name__ == "__main__":
    main()
  
    

