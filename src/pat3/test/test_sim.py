#! /usr/bin/env python
import os, math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.guidance as p3_guid
import pat3.atmosphere as p3_atm
import pat3.utils as p3_u

def main():
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)

    param_filename='/home/poine/work/pat/data/vehicles/cularis.xml'
    trim_args = {'h':0, 'va':10, 'gamma':0, 'flaps':0}
    _fdm = p1_fw_dyn.DynamicModel(param_filename)
    _fms = p3_guid.FMS(_fdm, trim_args)
    _atm = p3_atm.AtmosphereCalm()
    #_atm = p3_atm.AtmosphereThermalMulti()
    sim = p3_u.Sim(_fdm, _fms, _atm)
    
    t0, tf, dt = 0, 5., 0.01
    time = np.arange(t0, tf, dt)
    Xe, Ue = _fdm.trim(trim_args, report=True)

    sim.reset(time[0], Xe, Ue)
    sim.ctl.set_mode(p3_guid.FMS.mod_auto1)
    for i in range(1, len(time)):
        sim.run(time[i])
    _fdm.plot_trajectory(time, np.array(sim.Xs), np.array(sim.Us), window_title="None")
    plt.show()
    
if __name__ == "__main__":
    main()
