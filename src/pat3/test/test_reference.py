#!/usr/bin/env python3
import sys, time, numpy as np
import scipy.signal

import matplotlib.pyplot as plt

import pat3.utils as p3_u
import pat3.plot_utils as p3_pu

'''
Unit test for saturated linear reference model
'''

import pdb

def plot_ref(_time, Xr, _id="", _sp=None, fig=None, axs=None):
    #pdb.set_trace()
    fig, axs = plt.subplots(Xr.shape[1], 1) if fig is None else (fig, axs)
    for i in range(Xr.shape[1]):
        axs[i].plot(_time, Xr[:,i], label=_id)
        p3_pu.decorate(axs[i], f'$x^{{({i})}}$')
    if _sp is not None: axs[0].plot(_time, _sp, label='setpoint')
    p3_pu.decorate(axs[0], 'x', legend=True)
    return fig, axs

def run_ref(_ref, _time, _sp):
    Xr = np.zeros((len(_time), _ref.order+1))
    for i in range(1, len(_time)):
        Xr[i] = _ref.run( _time[i] - _time[i-1], _sp[i])
    return Xr

def _sec_order():
    poles, sats = p3_u.omxi_to_lambda(6., 0.9), [5., 25.]  # vel, accel
    return poles, sats

def _third_order():
    poles, sats = [-4, -4+2j, -4-2j], [6., 25., 50.] # vel accel jerk
    return poles, sats

def _fourth_order():
    poles = p3_u.omxi_to_lambda(3, 0.9)+p3_u.omxi_to_lambda(6, 0.8)
    sats = [6., 25., 50., 500.] # vel accel jerk jerkd
    return poles, sats

def _fifth_order():
    poles = p3_u.omxi_to_lambda(3, 0.9)+p3_u.omxi_to_lambda(6, 0.8)+ (-10.,)
    sats = [6., 25., 50., 500., 2500.] # vel accel jerk jerkd
    return poles, sats
    
def test(_ref_params, show_lin=False):
    poles, sats = _ref_params()
    K = -np.real(np.polynomial.polynomial.polyfromroots(poles)[:-1])
    _ref1 = p3_u.LinRef(K, sats)
    _time = np.arange(0, 10, 0.01)
    _sp = 4.*scipy.signal.square(_time*np.pi/3)
    Xr1 = run_ref(_ref1, _time, _sp)
    fig, axs = plot_ref(_time, Xr1, _id="saturated")
    if show_lin:
        _ref2 = p3_u.LinRef(K, None)
        Xr2 = run_ref(_ref2, _time, _sp)
        plot_ref(_time, Xr2, _id="linear", _sp=_sp, fig=fig, axs=axs)
    plt.show()
    
def main(args):
    test(_sec_order, True)
    #test(_third_order)
    #test(_fourth_order)
    #test(_fifth_order)
      
if __name__ == '__main__':
    main(sys.argv)
