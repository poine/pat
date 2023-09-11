#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

import pdb


def test1():
    x = np.array([0, 1, 0])
    l = np.array([0, 1, 2])

    l = np.array([ 0. ,  1.2,  1.9,  4.2,  3. ,  1.5])
    x = np.array([ 0. ,  0.2,  0.4 ,  0.6,  0.8,  1.])
    
    t, c, k = interpolate.splrep(x, l, s=0, k=4)
    print(f'''\
    t: {t}
    c: {c}
    k: {k}
    ''')
    N = 100
    xmin, xmax = x.min(), x.max()
    xx = np.linspace(xmin, xmax, N)
    spline = interpolate.BSpline(t, c, k, extrapolate=False)
    
    plt.plot(x, l, 'bo', label='Original points')
    plt.plot(xx, spline(xx), 'r', label='BSpline')


def test2():
    x = np.array([0, 1, 0, 1, 0, 1, 1.5])
    y = np.array([0, 1, 2, 3, 4, 5, 0])
    l = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.])

    sx = interpolate.InterpolatedUnivariateSpline(l, x, k=4)
    sy = interpolate.InterpolatedUnivariateSpline(l, y, k=4)
    l1 = np.arange(0, 1., 0.01)
    x1 = sx(l1)
    y1 = sy(l1)


    
    nder = 5
    Yc = np.zeros((2, nder, len(l1)))
    Yc[0] = np.array([sx.derivatives(_l) for _l in l1]).T
    Yc[1] = np.array([sy.derivatives(_l) for _l in l1]).T
       
    
    fig, axs = plt.subplots(nder, 2)
    axs[0,0].plot(l, x, 'bo')
    #axs[0,0].plot(l1, x1)
    axs[0,1].plot(l, y, 'bo')
    #axs[0,1].plot(l1, y1)
    for i in range(nder):
        axs[i,0].plot(l1, Yc[0, i])
        axs[i,1].plot(l1, Yc[1, i])

    
    plt.figure()
    plt.plot(x, y, 'bo')
    plt.plot(x1, y1)

    

test2()
plt.show()
