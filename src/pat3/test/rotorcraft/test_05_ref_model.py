#!/usr/bin/env python3
import sys, time, numpy as np
import scipy.signal
import matplotlib.pyplot as plt

import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt

import pdb


def main(args):
    #traj = pmt.RefModTraj([0, 0, -0.25], [1, 1, -0.25], v=0.1)
    if 0: # second order
        nder=3
        omega, xi = 6., 0.9
        poles = p3_u.omxi_to_lambda(omega, xi)
        print(poles)
        K = -np.real(np.polynomial.polynomial.polyfromroots(poles)[:-1])
        print(K)
        #K = np.array([-omega**2, -2*xi*omega])
        #print("good", K)
        sats = [6., 25.] # vel accel
        #sats = None
    if 0: # 3rd order
        nder = 4
        poles = [-4, -4+2j, -4-2j]
        K = -np.real(np.polynomial.polynomial.polyfromroots(poles)[:-1])
        sats = [6., 25., 50.] # vel accel jerk
        #sats = None
    if 1:
        nder=5
        poles = p3_u.omxi_to_lambda(3, 0.9)+p3_u.omxi_to_lambda(6, 0.8)
        #poles = [-1, -2, -3, -4]
        K = -np.real(np.polynomial.polynomial.polyfromroots(poles)[:-1])
        sats = [6., 25., 50., 500.] # vel accel jerk jerkd
        #sats = None
        print(poles)
        print(K)
  
    _ref = p3_u.LinRef(K, sats)
    _time = np.arange(0, 10, 0.01) 
    _sp = 4.*scipy.signal.square(_time*np.pi/3)
    Xr = np.zeros((len(_time), nder))
    for i in range(1, len(_time)):
        dt = _time[i] - _time[i-1]
        Xr[i] = _ref.run(dt, _sp[i])
    
    #pdb.set_trace()
    fig, axs = plt.subplots(nder, 1)
    for i in range(nder):
        axs[i].plot(_time, Xr[:,i])
    axs[0].plot(_time, _sp)
    plt.show()
        
if __name__ == '__main__':
    main(sys.argv)
