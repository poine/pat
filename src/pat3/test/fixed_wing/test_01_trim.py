#! /usr/bin/env python
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import control

''' Testing fixed wing trims '''

import pat3.utils as p3_u
import pat3.frames as p3_fr
import pat3.dynamics as pat_dyn
import pat3.vehicles.fixed_wing.simple_6dof_fdm as fw_dyn
import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn

def get_default_dm(ac_name='cularis'):#ac_name='skywalker_x8'):
    return p1_fw_dyn.DynamicModel(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))
    #return p1_fw_dyn.DynamicModel_ee(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))

def trimed_vz():
    dm = get_default_dm()
    v, R, g = 9., 15., 9.81
    Xe, Ue = dm.trim({'h':0., 'va':v, 'throttle':0}, debug=True)
    for R in [15. ]:
        #v^2 = R.g.tan(phi)
        phi = np.arctan(v**2/g/R)
        print('v:{}m/s R:{}m phi:{:.1f}deg'.format(v, R, np.rad2deg(phi)))
        gamma_e, ele_e, alpha_e = p1_fw_dyn.trim_banked(dm.P,  h=0., va=v, throttle=0., phi=phi, report=True, debug=True)
        #gamma_e, ele_e, alpha_e = p1_fw_dyn.trim_cst_throttle(dm.P,  h=0., va=v, throttle=0., report=True, debug=True)
        #print gamma_e, ele_e, alpha_e


        

def get_trim_defaults(trim_args={'h':0, 'va':12, 'gamma':0}, ac_name='cularis'):
    dm = p1_fw_dyn.DynamicModel(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))
    Xe, Ue = dm.trim(trim_args, debug=True)
    return dm, Xe, Ue

def print_trim(trim_args={'h':0, 'va':12, 'gamma':0}):
    dm, Xe, Ue = get_trim_defaults(trim_args)
    print(dm.state_str(Xe))

def format_evec(_evec, eps=np.finfo(float).eps, tol=100, prec=1):
    txt = '\n'
    for _icmp, _cmp in enumerate(_evec):
        txt += '{} '.format(np.array2string(np.real(_cmp) if np.imag(_cmp)<tol*eps else _cmp, precision=prec))
        if _icmp%3==2: txt += "    "
    return txt
    
def print_trim_dyn(trim_args={'h':0, 'va':12, 'gamma':0}, ac_name='cularis'):
    dm, Xe, Ue = get_trim_defaults(trim_args, ac_name)
    A, B = dm.get_jacobian(Xe, Ue)
    _evals, _evecs = np.linalg.eig(A)
    print('Modes:')
    for i, (_eval, _evec) in enumerate(zip(_evals, _evecs)):
        if np.isreal(_eval): # real eigenvalue
            txt = 'lambda_{}  : {:.3f}'.format(i+1, np.real(_eval))
            if np.abs(np.real(_eval))>1e-6: txt += ' (tau {:.1f}s)'.format(np.abs(1./np.real(_eval)))
            txt += format_evec(_evec)
            print(txt)
        else:               # complex eigenvalue
            if i>0 and np.imag(_eval)==-np.imag(_evals[i-1]):
                txt = ''    # conjugate has already been printed
            else:
                om = np.linalg.norm(_eval); xi = -np.real(_eval)/om
                txt = 'lambda_{}/{}: {:.3f} (om {:.2f} rad/s, xi {:.2f})'.format(i+1, i+2, _eval, om, xi)
                txt += format_evec(_evec)
                print(txt)
            

def plot_poles(trim_args={'h':0, 'va':12, 'gamma':0}, ac_name='cularis'):
    dm = p1_fw_dyn.DynamicModel(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))
    Xe, Ue = dm.trim(trim_args, debug=True)
    A, B = dm.get_jacobian(Xe, Ue)
    _eval, _evel = np.linalg.eig(A)
    B1 = B[:,2][np.newaxis].T
    C = np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0])
    D = np.zeros((1,1))
    ss = control.ss(A, B1, C, D)
    tf = control.ss2tf(ss)
    pdb.set_trace()


def test01(trim_args={'h':0, 'va':15, 'throttle':0, 'flaps':np.deg2rad(-8.)}, ac_name='cularis'):
    dm = p1_fw_dyn.DynamicModel(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))
    #dm = p1_fw_dyn.DynamicModel_ee(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))
    #dm = p1_fw_dyn.DynamicModel_eq(os.path.join(p3_u.pat_dir(), 'data/vehicles/{}.xml'.format(ac_name)))
    Xe, Ue = dm.trim(trim_args, debug=True, report=True, aero=False)
    print Xe, Ue
    Xe, Ue = dm.trim(trim_args, debug=True, report=True, aero=True)
    print Xe, Ue
    return dm, Xe, Ue
    
def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    test01()
    #test_equilibrium() # for now-broken fdm
    #trimed_vz() # trimming in circle...
    #print_trim()
    #print_trim_dyn(trim_args={'h':0, 'va':20, 'gamma':0}, ac_name='skywalker_x8')
    #plot_poles()
    plt.show()


    
if __name__ == "__main__":
    main()
