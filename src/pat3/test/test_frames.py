#! /usr/bin/env python
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import pat3.frames as p3_fr

def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)

    n_eulers = 1000
    eulers = np.random.uniform([-np.pi, -np.pi/2, -np.pi], [np.pi, np.pi/2, np.pi], (n_eulers, 3))
    eulers = np.array([[0, 0, 0], [0, 0, np.deg2rad(45)]])
    for euler in eulers:
        pos_ned = [0, 0, 0]
        ivel_ned = [10, 0, 0]
        avel_aero = va, alpha, beta = p3_fr.vel_world_to_aero(pos_ned, ivel_ned, euler, atm=None)
        ivel_ned2 = p3_fr.vel_aero_to_world(pos_ned, avel_aero, euler, atm=None)
        pdb.set_trace()
        
        print np.allclose(ivel_ned, ivel_ned2)
    

if __name__ == "__main__":
    main()
