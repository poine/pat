#!/usr/bin/env python3
import sys, time, numpy as np
import scipy.signal
#import rospy
import matplotlib.pyplot as plt

import pat3.utils as p3_u
import pat3.plot_utils as p3_pu

def plot_ref(_time, Xr, _id="", _sp=None, fig=None, axs=None):
    fig, axs = plt.subplots(3, 1) if fig is None else (fig, axs)
    axs[0].plot(_time, Xr[:,0], label=_id)
    if _sp is not None: axs[0].plot(_time, _sp, label='setpoint')
    p3_pu.decorate(axs[0], 'pos', legend=True)
    axs[1].plot(_time, Xr[:,1])
    p3_pu.decorate(axs[1], 'vel')
    axs[2].plot(_time, Xr[:,2])
    p3_pu.decorate(axs[2], 'accel')
    return fig, axs

def run_ref(_ref, _time, _sp):
    Xr = np.zeros((len(_time), 3))
    for i in range(1, len(_time)):
        dt = _time[i] - _time[i-1]
        Xr[i] = _ref.run(dt, _sp[i])
    return Xr
    
def test_scalar_ref():
    _time = np.arange(0, 10, 0.01)

    #_sats = [6., 50.]  # vel, accel
    _sats = [6., 25.]  # vel, accel
    _ref1 = p3_u.SecOrdLinRef(omega=6, xi=0.9, sats=None)
    _ref2 = p3_u.SecOrdLinRef(omega=6, xi=0.9, sats=_sats)

    #_sp = 1.*np.ones(len(_time))
    _sp = 4.*scipy.signal.square(_time*np.pi/3)
    Xr1 = run_ref(_ref1, _time, _sp)
    Xr2 = run_ref(_ref2, _time, _sp)

    fig, axs = plot_ref(_time, Xr1, _id="linear")
    plot_ref(_time, Xr2, _id="saturated", _sp=_sp, fig=fig, axs=axs)
    plt.show()
        
def main(args):
    test_scalar_ref()
    
      
if __name__ == '__main__':
    main(sys.argv)
