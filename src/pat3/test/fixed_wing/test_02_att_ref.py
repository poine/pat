#! /usr/bin/env python
import os, sys, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb


'''
  Testing attitude reference model
'''

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil
import pat3.utils as p3_u
import pat3.atmosphere as p3_atm
import pat3.frames as p3_fr
import pat3.plot_utils as p3_pu



def main(t0=0, tf=5., dt=0.005):
    va = 12.
    time = np.arange(t0, tf, dt)
    ref = p3_pil.EulerRef()
    ref_logger = p3_pil.EulerRefLogger()
    phi_sp, theta_sp = np.zeros(len(time)), np.zeros(len(time))
    theta_sp = np.deg2rad(2)+np.array([p3_u.step(_t, a=np.deg2rad(10), p=10., dt=2.5) for _t in time])
    phi0, theta0 = np.deg2rad(10), 0
    ref.reset(phi0, theta0)
    ref_logger.log(ref)
    for i in range(1, len(time)):
        ref.run(dt, phi_sp[i], theta_sp[i], va)
        ref_logger.log(ref)
    ref_logger.plot(time)
    plt.show()
    
    
        
   
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    main()
