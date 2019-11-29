#! /usr/bin/env python
#-*- coding: utf-8 -*-
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

''' Checking match coefficients between AVL and XFLR5 '''
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.frames as p3_fr
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.simple_6dof_fdm_param


def XFLR5_to_AVL(P):
    # overwrite (unknown) AVL coefficients with values computed from XFLR5
    print('old pat: CD0 {} CD_k1 {} CD_k2 {}'.format(P.CD0, P.CD_k1, P.CD_k2))
    CD0, CD_a1, CD_a2 = 0.0107, -0.00955, 1.1  # XFLR5 coefs
    P.CD_k2 = CD_a2/P.CL_alpha**2              # AVL counterparts
    P.CD_k1 = CD_a1/P.CL_alpha - 2*P.CD_k2*P.CL0
    P.CD0 = CD0 - P.CD_k1*P.CL0 - P.CD_k2*P.CL0**2
    print('new pat: CD0 {} CD_k1 {} CD_k2 {}'.format(P.CD0, P.CD_k1, P.CD_k2))
    
def get_x8_aero_coef(alpha, beta, rvel, Usfc, P):
    d_alpha = alpha - P.alpha0
    nrvel = rvel*np.array([P.Bref, P.Cref, P.Bref])/2/P.Vref
    de = Usfc[1]
    CL = P.CL0 + P.CL_alpha*d_alpha + np.dot(P.CL_omega,nrvel) + np.dot(P.CL_sfc,Usfc)
    CY = 0
    CD0, CDa1, CDa2 = 0.0107, -0.00955, 1.1
    CDde2, CDbeta2 = 0.0196, 0.115
    CD = CD0 + CDa1*d_alpha + CDa2*d_alpha**2 + CDde2*de**2 + CDbeta2*beta**2
    return [CL, CY, CD]


def get_both_coefs(P, alpha, beta, rvel, Usfc):
    CL1, CY1, CD1 = p1_fw_dyn.get_f_aero_coef(alpha, beta, rvel, Usfc, P)
    CL2, CY2, CD2 = get_x8_aero_coef(alpha, beta, rvel, Usfc, P)
    return (CL1, CY1, CD1), (CL2, CY2, CD2)

def test_alpha(P):
    # alpha variation (de, rvel = 0)
    alpha, beta, rvel = 0., 0., [0., 0., 0.]
    Usfc = [0, 0, 0, 0]
    alphas = np.deg2rad(np.linspace(-2, 15, 100))
    coefs = np.array([get_both_coefs(P, alpha, beta, rvel, Usfc) for alpha in alphas])
    ax = plt.subplot(3,2,1)
    plt.plot(np.rad2deg(alphas), coefs[:,0,0], label='CL avl')
    plt.plot(np.rad2deg(alphas), coefs[:,1,0], label='CL XFLR5')
    p3_pu.decorate(ax, title='CL', xlab='alpha (deg)', ylab=None, legend=True)
    ax = plt.subplot(3,2,2)
    plt.plot(np.rad2deg(alphas), coefs[:,0,2], label='CD avl')
    plt.plot(np.rad2deg(alphas), coefs[:,1,2], label='CD XFLR5')
    p3_pu.decorate(ax, title='CD', xlab='alpha (deg)', ylab=None, legend=True)
    # elevator variation (alpha, rvel = 0)
    deles =  np.deg2rad(np.linspace(-15, 15, 100))
    Usfcs = np.zeros((len(deles), 4)); Usfcs[:,1] = deles
    coefs = np.array([get_both_coefs(P, alpha, beta, rvel, Usfc) for Usfc in Usfcs])
    ax = plt.subplot(3,2,3)
    plt.plot(np.rad2deg(deles), coefs[:,0,0], label='CL avl')
    plt.plot(np.rad2deg(deles), coefs[:,1,0], label='CL XFLR5')
    p3_pu.decorate(ax, title='CL', xlab='d_ele (deg)', ylab=None, legend=True)
    ax = plt.subplot(3,2,4)
    plt.plot(np.rad2deg(deles), coefs[:,0,2], label='CD avl')
    plt.plot(np.rad2deg(deles), coefs[:,1,2], label='CD XFLR5')
    p3_pu.decorate(ax, title='CD', xlab='d_ele (deg)', ylab=None, legend=True)
    # pdb.set_trace()



    plt.show()

def main():
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/skywalker_x8.xml')
    P = pat3.vehicles.fixed_wing.simple_6dof_fdm_param.Param(param_filename)
    XFLR5_to_AVL(P)
    print(P)
    test_alpha(P)

    
    alpha, beta, rvel = 0., 0., [0., 0., 0.]
    Usfc = [0, 0, 0, 0]
    (CL1, CY1, CD1), (CL2, CY2, CD2) = get_both_coefs(P, alpha, beta, rvel, Usfc)
    print(CL1, CY1, CD1)
    print(CL2, CY2, CD2)
    
if __name__ == "__main__":
    main()

