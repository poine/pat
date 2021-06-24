#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os, os.path, pickle, logging, numpy as np, matplotlib.pyplot as plt
import scipy.interpolate
#import keras
import control
import pdb

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.utils as p3_u
import pat3.plot_utils as p3_pu
import pat3.frames as p3_fr

'''

  Playing with gain scheduling

'''
#
# TODO: finish that :)
# 
def _get_longi_lti_ssr(dm, trim_args, dt=0.01, report=False):
    Xe, Ue = dm.trim(trim_args, report=report)
    Ac, Bc = dm.get_jacobian(Xe, Ue)
    s = p3_fr.SixDOFAeroEuler
    A1c = np.array([[Ac[s.sv_theta, s.sv_theta], Ac[s.sv_theta, s.sv_q]],
                    [Ac[s.sv_q, s.sv_theta], Ac[s.sv_q, s.sv_q]]])
    B1c = np.array([[Bc[s.sv_theta, dm.iv_de()]],
                    [Bc[s.sv_q,     dm.iv_de()]]])
    ct_ssr = control.ss(A1c, B1c, [[1, 0]], [[0]])
    dt_ssr = control.sample_system(ct_ssr, dt, method='zoh') #  ‘matched’, ‘tustin’, ‘zoh’
    return A1c, B1c, dt_ssr.A, dt_ssr.B
    

class ScheduledPitchControl:
    def __init__(self, dm, om_eps=10., xi_eps=0.7, om_r=5., xi_r=0.9, dt=0.01):
        self.ref = p3_u.SecOrdLinRef(omega=om_r, xi=xi_r, sats=[np.deg2rad(45.), np.deg2rad(100.)])  # theta, q
        self.compute_gains(dm, om_eps, xi_eps, om_r, xi_r, dt)

    def compute_gains(self, dm, om_eps, xi_eps, om_r, xi_r, dt):
        vas = np.arange(5., 30., 1.)
        a0s, a1s, bs = [], [], []
        for va in vas:
            trim_args = {'h':0, 'va':va, 'gamma':0.}
            Ac, Bc, Ad, Bd = _get_longi_lti_ssr(dm, trim_args, dt)
            a0s.append(-Ac[1,0]); a1s.append(-Ac[1,1]); bs.append(Bc[1,0])
            # we could use place
            #Ks.append(control.place(Ac, Bc, cl_poles))
        
        om2_eps, txiom_eps = om_eps**2, 2.*xi_eps*om_eps 
        Ks = np.array([[-(a0-om2_eps)/b, -(a1-txiom_eps)/b] for a0, a1, b in zip(a0s, a1s, bs)])
        om2_r, txiom_r = om_r**2, 2.*xi_r*om_r 
        Hs = np.array([[(om2_eps-om2_r)/b, (txiom_eps-txiom_r)/b, om2_r/b] for b in bs])

        self.k0_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, Ks[:,0])
        self.k1_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, Ks[:,1])
        self.h0_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, Hs[:,0])
        self.h1_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, Hs[:,1])
        self.h2_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, Hs[:,2])

    def _gains(self, va):
        return [self.k0_spl(va), self.k1_spl(va)], [self.h0_spl(va), self.h1_spl(va), self.h2_spl(va)]

    def get(self, X, t, theta_sp):
        #theta_refkp1, q_refkp1, qd_refkp1 = self.ref.run(self.dt, theta_sp-self.theta_e)
        return 0.

    

def run_sim(plant, ctl):
    trim_args = {'h':0, 'va':12, 'gamma':0.}
    Xe, Ue = plant.trim(trim_args, report=True)
    _s = p3_fr.SixDOFAeroEuler
    time =  np.arange(0., 30.01, plant.dt)
    X, U = np.zeros((len(time), _s.sv_size)), np.zeros((len(time), plant.input_nb()))
    sp = p3_u.step_vec(time, a=np.deg2rad(-1), dt=20)
    X[0] = np.array(Xe); X[0, _s.sv_theta] += np.deg2rad(2)
    plant.reset(X[0])
    for k in range(0,len(time)-1):
        U[k] = ctl.get(X[k], time[k], sp[k])
        X[k+1] = plant.run(plant.dt, time[k+1], U[k], atm=None)
        
        
def main():
    param_filename = os.path.join(p3_u.pat_dir(), 'data/vehicles/cularis.xml')
    dm = p1_fw_dyn.DynamicModel(param_filename)

    pc = ScheduledPitchControl(dm)
    run_sim(dm, pc)

    #old_thing(dm)
    plt.show()

def old_thing(dm):
    om_eps, xi_eps = 10., 0.7  # perturbation rejection
    om_r, xi_r = 5., 0.9       # reference model
    cl_poles = p3_u.omxi_to_lambda(om_eps, xi_eps)
    #vas = [5., 10., 15., 20., 30.]
    vas = np.arange(5., 30., 1.)
    a0s, a1s, bs = [], [], []
    Ks = []
    for va in vas:
        trim_args = {'h':0, 'va':va, 'gamma':0.}
        Ac, Bc, Ad, Bd = _get_longi_lti_ssr(dm, trim_args)
        a0s.append(-Ac[1,0]); a1s.append(-Ac[1,1]); bs.append(Bc[1,0]) 
        Ks.append(control.place(Ac, Bc, cl_poles))
    Ks = np.array(Ks).squeeze()
    
    om2_eps, txiom_eps = om_eps**2, 2.*xi_eps*om_eps 
    K2s = np.array([[-(a0-om2_eps)/b, -(a1-txiom_eps)/b] for a0, a1, b in zip(a0s, a1s, bs)])
    om2_r, txiom_r = om_r**2, 2.*xi_r*om_r 
    Hs = np.array([[(om2_eps-om2_r)/b, (txiom_eps-txiom_r)/b, om2_r/b] for b in bs])


    #pdb.set_trace()
    ax = plt.subplot(3,3,1)
    plt.plot(vas, a0s)
    p3_pu.decorate(ax, title='a0', xlab='va (m/s)')
    ax = plt.subplot(3,3,2)
    plt.plot(vas, a1s)
    p3_pu.decorate(ax, title='a1', xlab='va (m/s)')
    ax = plt.subplot(3,3,3)
    plt.plot(vas, bs)
    p3_pu.decorate(ax, title='b', xlab='va (m/s)')

    ax = plt.subplot(3,3,4)
    plt.plot(vas, Ks[:,0])
    plt.plot(vas, K2s[:,0])
    p3_pu.decorate(ax, title='k0', xlab='va (m/s)')
    ax = plt.subplot(3,3,5)
    plt.plot(vas, Ks[:,1])
    plt.plot(vas, K2s[:,1])
    p3_pu.decorate(ax, title='k1', xlab='va (m/s)')

    ax = plt.subplot(3,3,7)
    plt.plot(vas, Hs[:,0])
    p3_pu.decorate(ax, title='h0', xlab='va (m/s)')
    ax = plt.subplot(3,3,8)
    plt.plot(vas, Hs[:,1])
    p3_pu.decorate(ax, title='h1', xlab='va (m/s)')
    ax = plt.subplot(3,3,9)
    plt.plot(vas, Hs[:,2])
    p3_pu.decorate(ax, title='h2', xlab='va (m/s)')


    # plt.subplot(3,3,6)
    # plt.plot(vas, K2s[:,0]*vas**2)

    if 1:
        #spl = scipy.interpolate.UnivariateSpline(vas, K2s[:,0])
        k0_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, K2s[:,0])
        k1_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, K2s[:,1])
        h0_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, Hs[:,0])
        h1_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, Hs[:,1])
        h2_spl = scipy.interpolate.InterpolatedUnivariateSpline(vas, Hs[:,2])
        plt.subplot(3,3,4)
        plt.plot(vas, k0_spl(vas))
        plt.subplot(3,3,5)
        plt.plot(vas, k1_spl(vas))
        plt.subplot(3,3,7)
        plt.plot(vas, h0_spl(vas))

    

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    main()
