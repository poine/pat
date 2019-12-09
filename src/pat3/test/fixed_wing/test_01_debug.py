#! /usr/bin/env python
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

''' Testing fixed wing in open loop '''
import pat3.utils as p3_u
import pat3.frames as p3_fr
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.atmosphere as p3_atm

import test_01_dynamics as t1dyn


def test1():
    atm = p3_atm.AtmosphereCstWind([0.5, 0, 0.])
    dm = t1dyn.get_default_dm()
    time = np.arange(0, 30, 0.01)
    Xe, Ue = dm.trim({'h':0, 'va':10, 'throttle':0}, debug=True)
    X0, X_act0 = np.array(Xe), np.array(Ue),
    U = Ue*np.ones((len(time), dm.input_nb()))
    time, X, X_act = t1dyn.run_simulation(dm, time, X0, X_act0, U, atm=atm)
    dm.plot_trajectory(time, X, X_act, window_title='aero/euler')
    dm.plot_trajectory_ee(time, X, X_act, window_title='euclidian/euler', atm=atm)
    plt.show()


def test_ee_vs_ae():
    atm = None#p3_atm.AtmosphereCstWind([0.5, 0., 0.])
    dm_ee = p1_fw_dyn.DynamicModel_ee()
    dm_ae = p1_fw_dyn.DynamicModel()
    if 1:
        U = [0., 0.0, 0, 0, 0]
        Xee = np.array([0., 0., 0.,  10., 0., 0.,  0., 0., 0.,  0., 0., 0.])
        print('Xee {}'.format(Xee))
        Xee_dot = dm_ee.dyn(Xee, t=0, U=U, atm=atm)
        print('Xee_dot {}'.format(Xee_dot))
        Xae = p3_fr.SixDOFEuclidianEuler.to_six_dof_aero_euler(Xee, atm)
        print('Xae {}'.format(Xae))
        Xae_dot = dm_ae.dyn(Xae, t=0, U=U, atm=atm)
        print('Xae_dot {}'.format(Xae_dot))
    if 0:
        #Xae_e, Uae_e = dm_ae.trim({'h':0, 'va':10, 'throttle':0}, debug=True)
        Xae_e, Uae_e = dm_ae.trim({'h':0, 'va':10, 'gamma':0}, report=True, debug=True)
        print('Xae_e {}'.format(Xae_e))
        print
        #Xee_e, Uee_e = dm_ee.trim({'h':0, 'va':10, 'throttle':0}, debug=True)
        Xee_e, Uee_e = dm_ee.trim({'h':0, 'va':10, 'gamma':0}, report=True, debug=True)
        print('Xee_e {}'.format(Xee_e))

        time = np.arange(0, 10, 0.01)

        X0, X_act0 = np.array(Xee_e), np.array(Uee_e),
        U = Uee_e*np.ones((len(time), dm_ee.input_nb()))
        time, Xee, Xee_act = t1dyn.run_simulation(dm_ee, time, X0, X_act0, U, atm=atm)
        dm_ee.plot_trajectory_ee(time, Xee, Xee_act, window_title='euclidian/euler')
    
        X0, X_act0 = np.array(Xae_e), np.array(Uae_e),
        U = Uee_e*np.ones((len(time), dm_ae.input_nb()))
        time, Xae, Xae_act = t1dyn.run_simulation(dm_ae, time, X0, X_act0, U, atm=atm)
        dm_ae.plot_trajectory_ee(time, Xae, Xae_act, window_title='aero/euler')
    
        plt.show()
    
def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    #test1()
    test_ee_vs_ae()

if __name__ == "__main__":
    main()
