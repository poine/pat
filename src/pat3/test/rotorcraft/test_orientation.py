#! /usr/bin/env python



import numpy as np
import pdb

import pat3.transformations as ptr
import pat3.algebra as pal





def test_quaternion_from_euler():
    nb_samples = 1000
    phis = np.random.uniform(low=-np.pi, high=np.pi, size=nb_samples)
    thetas = np.random.uniform(low=-np.pi/2, high=np.pi/2, size=nb_samples)
    psis = np.random.uniform(low=-np.pi, high=np.pi, size=nb_samples)

    n_fail = 0
    for phi, theta, psi in zip(phis, thetas, psis):
        q_tr  = ptr.quaternion_from_euler(psi, theta, phi, 'rzyx')
        q_pal = pal.quat_of_euler([phi, theta, psi])
        if not np.allclose(q_tr, pal.q_ixyz_to_xyzw(q_pal)):
            n_fail += 1
    print('quaternion_from_euler {} failure ({} samples)'.format(n_fail, nb_samples))


def main():
    test_quaternion_from_euler()


    

if __name__ == "__main__":
    np.set_printoptions(linewidth=500)
    main()
