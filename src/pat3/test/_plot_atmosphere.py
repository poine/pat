#!/usr/bin/env python
#-*- coding: utf-8 -*-
#

import numpy as np, matplotlib.pyplot as plt
import matplotlib, matplotlib.pyplot as plt, mpl_toolkits.mplot3d
from matplotlib import cm


import pat3.plot_utils as p3_plu, pat3.atmosphere as p3_atm

def display_shear():
    atm1 =  p3_atm.AtmosphereShearX(wind1=5.0, wind2=0.0, xlayer=60.0, zlayer=40.0)
    atm2 =  p3_atm.AtmosphereVgradient(w0=0, w1=5, h0=20, h1=60)
    #atm2 =  p3_atm.AtmosphereVgradient(w0=-7.5, w1=7.5, h0=20, h1=60)
    figure = plt.figure(figsize=(10.24, 3.12))
    axes = figure.subplots(1, 2)
    p3_plu.plot_slice_wind_nu(atm1, n0=-100, n1=80, dn=10., e0=0., h0=-10, h1=100, dh=5., zdir=-1.,
                              show_quiver=True, show_color_bar=False,
                              title='Shear',
                              figure=figure, ax=axes[0], use_wx=True)
    p3_plu.plot_slice_wind_nu(atm2, n0=-100, n1=80, dn=10., e0=0., h0=-10, h1=100, dh=5., zdir=-1.,
                              show_quiver=True, show_color_bar=False,
                              title='Vgrad',
                              figure=figure, ax=axes[1], use_wx=True)
    plt.savefig('/tmp/atm_shear_vs_vgrad.png', dpi=120, bbox_inches='tight')


def display_all_models(h0=0.):
    atm1 =  p3_atm.AtmosphereWharington()#cst=[1.5, 0, 0])
    #atm2 =  p3_atm.AtmosphereShearX()
    atm2 = p3_atm.AtmosphereRidge()
    nc_f = '/home/poine/work/glider_experiments/data/extr_IHODC.1.RK4DI.007.nc'
    atm3 = p3_atm.AtmosphereNC(nc_f)
    atms, names = [atm1, atm2, atm3], ['Thermal', 'Shear', 'Grid']
    matplotlib.rcParams['text.usetex'] = True
    plt.rcParams["font.family"] = "Times New Roman"  
    plt.rcParams["font.size"] = 11  
    #fig = plt.figure(figsize=(10,5))
    figure = plt.figure(figsize=(10.24, 3.12))
    axes = figure.subplots(1, 3)
    #plt.tight_layout()
    left, bottom, right, top, wspace, hspace = 0.07, 0.14, 0.97, 0.89, 0.26, 0.3
    figure.subplots_adjust(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)
    #
    p3_plu.plot_slice_wind_nu(atm1, n0=-80, n1=80, dn=10., e0=0., h0=-10, h1=100, dh=5., zdir=-1.,
                              show_quiver=True, show_color_bar=False,
                              title='Wharington',
                              figure=figure, ax=axes[0])
    #
    p3_plu.plot_slice_wind_nu(atm2, n0=-60, n1=100, dn=10., e0=0., h0=-10, h1=100, dh=5., zdir=-1.,
                              show_quiver=True, show_color_bar=False,
                              title='Ridge',
                              figure=figure, ax=axes[1])
    #
    p3_plu.plot_slice_wind_ne(atm3, n0=0, n1=1000, dn=50., e0=0., e1=1000, de=50, h0=100, t0=0.,
                              #show_quiver=True,
                              show_color_bar=False,
                              title='Grid',
                              figure=figure, ax=axes[2])

    plt.savefig('/tmp/atmosphere_examples.png')
 
def main():
    display_all_models()
    display_shear()
    plt.show()

    
if __name__ == "__main__":
    main()
  
    

