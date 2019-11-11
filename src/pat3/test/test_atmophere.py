#!/usr/bin/env python
#-*- coding: utf-8 -*-
#

import numpy as np, matplotlib.pyplot as plt
import matplotlib, matplotlib.pyplot as plt, mpl_toolkits.mplot3d
from matplotlib import cm
# This is for the thermal Bell shapes:

def thermal_model_allen(x,y,z,wstar=256,wgain=1.,rgain=1.,zi=1000):
    r1r2shape = [0.1400, 0.2500, 0.3600, 0.4700, 0.5800, 0.6900, 0.8000]
    Kshape = [[1.5352, 2.5826, -0.0113, -0.1950, 0.0008],
    [1.5265, 3.6054, -0.0176, -0.1265, 0.0005],
    [1.4866, 4.8356, -0.0320, -0.0818, 0.0001],
    [1.2042, 7.7904,  0.0848, -0.0445, 0.0001],
    [0.8816, 13.9720, 0.3404, -0.0216, 0.0001],
    [0.7067, 23.9940, 0.5689, -0.0099, 0.0002],
    [0.6189, 42.7965, 0.7157, -0.0033, 0.0001]]

    #CALCULATE AVERAGE UPDRAFT SIZE
    zzi=z/zi; rbar=(.102*zzi**(1/3))*(1-(.25*zzi))*zi;

    #CALCULATE AVERAGE UPDRAFT STRENGTH
    wtbar=(zzi**(1/3))*(1-1.1*zzi)*wstar;

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

    
    r=1
    rr2=r/r2

    # From this point we may use the simple Gedeon Gaussian model:
    R_ave = (r1+r2)/2.
    w = thermal_model_gedeon(x, y, z, R_t=R_ave, W_max=wc)

    # Or Allen's real model
    if z<zi:
        if  r1r2<.5*(r1r2shape[0]+r1r2shape[1]) : i = 0
        elif r1r2<.5*(r1r2shape[1]+r1r2shape[2]) : i = 1
        elif r1r2<.5*(r1r2shape[2]+r1r2shape[3]) : i = 2
        elif r1r2<.5*(r1r2shape[3]+r1r2shape[4]) : i = 3
        elif r1r2<.5*(r1r2shape[4]+r1r2shape[5]) : i = 4
        elif r1r2<.5*(r1r2shape[5]+r1r2shape[6]) : i = 5
        else: i=6
        ka=Kshape[i][0]
        kb=Kshape[i][1]
        kc=Kshape[i][2]
        kd=Kshape[i][4] #check this again FIXME
        inn=rr2;
        #CALCULATE SMOOTH VERTICAL VELOCITY DISTRIBUTION
        ws = (1./(1+(ka*np.abs(inn+kc))**kb))+kd*inn;
    #     print(ws)

        # FIXmE ws[np.nonzero(ws<0)]=0; # no neg updrafts
    else:
        ws = 0.

    return w


def thermal_model_gedeon(x, y, z, R_t=10, W_max=3):
    r = np.sqrt(x**2+y**2)
    top = (r/R_t)**2
    w = W_max*np.exp(-top)*(1-top)
    return w


def Z_points(x,y,z,fun):
    X, Y = np.meshgrid(x, y)
    zs = np.array(fun(np.ravel(X), np.ravel(Y), z))
    Z = zs.reshape(X.shape)
    return X,Y,Z

def Z_points_2(x,y,z,fun):
    i=0; j=0;
    Z = np.zeros((len(x),len(y)))
    X, Y = np.meshgrid(x, y)
    for x_ in x:
        j=0
        for y_ in y:
            Z[i,j] = fun(x_,y_, z)
            j +=1
        i +=1
    return X,Y,Z


def main():
    x = y = np.arange(-60.0, 60.0, 1)
    z=50

    X,Y,Z = Z_points(x,y,z,thermal_model_gedeon)

    fig = plt.figure(figsize=(15,13))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z,  rstride=1, cstride=1, alpha=1.0);ax.plot_surface(X, Y, Z+7, alpha=1.0);
    # ax.quiver(X, Y, 12, 0, 0,Z, alpha=1.0);
    ax.set_xlabel('East [m]');ax.set_ylabel('North [m]');ax.set_zlabel('Up [m]')
    
    # x = y = np.arange(-60.0, 60.0, 1)
    
    fig = plt.figure(figsize=(15,13))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, rstride=2, cstride=2, alpha=0.8)
    cset = ax.contour(X, Y, Z, zdir='z', offset=0, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='x', offset=0, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, Z, zdir='y', offset=0, cmap=cm.coolwarm)


    fig = plt.figure(figsize=(15,13))
    ax = fig.add_subplot(111, projection='3d')
    X,Y,Z = Z_points(x,y,z,thermal_model_allen)
    ax.plot_surface(X, Y, Z+z, alpha=1.0)
    z=150
    X,Y,Z = Z_points(x,y,z,thermal_model_allen)
    ax.plot_surface(X, Y, Z+z, alpha=1.0);
    
    z=450
    X,Y,Z = Z_points(x,y,z,thermal_model_allen)
    ax.plot_surface(X, Y, Z+z, alpha=1.0);
    
    z=650
    X,Y,Z = Z_points(x,y,z,thermal_model_allen)
    ax.plot_surface(X, Y, Z+z, alpha=1.0);
    cset = ax.contour(X, Y, Z, zdir='x', offset=-65, cmap=cm.coolwarm)
    
    ax.set_xlabel('East [m]');ax.set_ylabel('North [m]');ax.set_zlabel('Up [m]')
    # ax.view_init(elev=0., azim=0)
    plt.show()
    


if __name__ == "__main__":
        main()



