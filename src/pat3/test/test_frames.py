#! /usr/bin/env python
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import pat3.frames as p3_fr


def test0():
    pos_ned = [0, 0, 0]
    #vel_aero = [10, np.deg2rad(4), 0]
    vel_aero = [10, 0., -0.6]
    eulers = np.deg2rad([0, 0, 45])
    #eulers = [1.77592771, 0.99355044, 3.13596614]
    #eulers = [ 1.23060692,  0.14490205, -2.96400653]
    rvel_body = [0, 0, 0]
    
    Xae = np.concatenate((pos_ned, vel_aero, eulers, rvel_body))
    Xee = p3_fr.SixDOFAeroEuler.to_six_dof_euclidian_euler(Xae, atm=None)
    print Xae
    print Xee
    va, alpha, beta = Xae[p3_fr.SixDOFAeroEuler.sv_slice_vaero]
    wind_ned = [0, 0, 0]
    ivel_ned =  p3_fr.vel_aero_to_world(vel_aero, eulers, wind_ned)
    vaero_2 = p3_fr.vel_world_to_aero_eul(ivel_ned, eulers, wind_ned)
    print vel_aero
    print ivel_ned
    print vaero_2
    print("allclose {}".format(np.allclose(vel_aero, vaero_2, rtol=1e-03, atol=1e-03, equal_nan=False)))
    #pdb.set_trace()

def test1():
    n_eulers = 1000
    eulers = np.random.uniform([-np.pi, -np.pi/2, -np.pi], [np.pi, np.pi/2, np.pi], (n_eulers, 3))
    #eulers = np.array([[0, 0, np.deg2rad(45)]])
    #eulers = [[1.77592771, 0.99355044, 3.13596614]]
    wind_ned = [0, 0, 0]
    pass_cnt = 0
    for euler in eulers:
        pos_ned = [0, 0, 0]
        ivel_ned = [10, 0, 0]
        avel_aero = va, alpha, beta = p3_fr.vel_world_to_aero_eul(ivel_ned, euler, wind_ned)
        ivel_ned2 = p3_fr.vel_aero_to_world(avel_aero, euler, wind_ned)

        if not np.allclose(ivel_ned, ivel_ned2, rtol=1e-03, atol=1e-03, equal_nan=False):
            print 'failed', ivel_ned, euler, ivel_ned2
        else: pass_cnt += 1

    print 'passed {}/{}'.format(pass_cnt, len(eulers))



    
def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    #test0()
    test1()
    

if __name__ == "__main__":
    main()
