#! /usr/bin/env python
import math, numpy as np, matplotlib.pyplot as plt
import logging
import pdb

import sys
#sys.path.append('/mnt/mint17/home/poine/work/python-aerospace-toolbox/src') # PAT1 now along PAT3
sys.path.append('/home/poine/work/pat/src') # PAT3
import pat.vehicles.fixed_wing.dynamic_model_python_basic as p1_fw_dyn
import pat.vehicles.fixed_wing.control_3d as p1_fw_ctl
import pat3.vehicles.fixed_wing.simple_6dof_fdm as p3_fw_dyn
import pat3.utils as p3_u


def compare_open_loop(dm1, Xe1, Ue1, dm3, Xe3, Ue3, tf=0.5, dt=0.01, Uf=None):
    time = np.arange(0, tf, dt)
    X1 = np.zeros((len(time), dm1.sv_size))
    U1 = np.zeros((len(time),  dm1.input_nb()))
    X3 = np.zeros((len(time), dm3.sv_size))
    U3 = np.zeros((len(time),  dm3.iv_size))
    X3[0] = Xe3
    U3[:] = Ue3
    X1[0] = dm1.reset(Xe1)
    U1[:] = Ue1
    print(p1_fw_dyn.dyn(X1[0], 0, U1[0], dm1.P))
    A1, B1 = dm1.get_jacobian(Xe1, Ue1)
    print(dm3.cont_dyn(X3[0], 0, U3[0]))
    #pdb.set_trace()
    for i in range(1, len(time)):
        if Uf is not None:
            U1[i-1] += Uf(time[i-1])
            U3[i-1] += Uf(time[i-1])
        X1[i] = dm1.run(time[i] - time[i-1], U1[i-1])
        #X3[i] = dm3.disc_dyn(X3[i-1], time[i-1], U3[i-1], time[i]-time[i-1])
    #dm3.plot(time, X3, U3)
    dm1.plot_trajectory(time, X1, U1)
    plt.show()

def compare_step_ele(dm1, Xe1, Ue1, dm3, Xe3, Ue3):
    def Uf(t): return [0, 0, p3_u.step(t, a=np.deg2rad(-0.5), p=5., dt=0.), 0, 0]
    return compare_open_loop(dm1, Xe1, Ue1, dm3, Xe3, Ue3, tf=10.5, dt=0.01, Uf=Uf)

def compare_step_ailerons(dm1, Xe1, Ue1, dm3, Xe3, Ue3):
    def Uf(t): return [0, p3_u.step(t, a=np.deg2rad(-1.), p=2., dt=1.), 0, 0, 0]
    return compare_open_loop(dm1, Xe1, Ue1, dm3, Xe3, Ue3, tf=10.5, dt=0.01, Uf=Uf)

def compare_step_rudder(dm1, Xe1, Ue1, dm3, Xe3, Ue3):
    def Uf(t): return [0, 0, 0, p3_u.step(t, a=np.deg2rad(-1.), p=2., dt=1.), 0]
    return compare_open_loop(dm1, Xe1, Ue1, dm3, Xe3, Ue3, tf=10.5, dt=0.01, Uf=Uf)


def compare_step_flaps(dm1, Xe1, Ue1, dm3, Xe3, Ue3):
    def Uf(t): return [0, 0, 0, 0, p3_u.step(t, a=np.deg2rad(-1.), p=2., dt=1.)]
    return compare_open_loop(dm1, Xe1, Ue1, dm3, Xe3, Ue3, tf=10.5, dt=0.01, Uf=Uf)

    
def compare_trim(param_filename='/home/poine/work/pat/data/vehicles/cularis.xml',
                 trim_args = {'h':0, 'va':15, 'gamma':0}):
    dm1 = p1_fw_dyn.DynamicModel(param_filename)
    Xe1, Ue1 = dm1.trim(trim_args, debug=True)
    dm3 = p3_fw_dyn.FDM3(p3_fw_dyn.Param(param_filename))
    Xe3, Ue3 = dm3.trim(trim_args, debug=True)
    return dm1, Xe1, Ue1, dm3, Xe3, Ue3 

def main():
    param_filename='/home/poine/work/pat/data/vehicles/cularis.xml'
    trim_args = {'h':0, 'va':15, 'gamma':0}
    dm1, Xe1, Ue1, dm3, Xe3, Ue3 = compare_trim(param_filename, trim_args)
    #compare_open_loop(dm1, Xe1, Ue1, dm3, Xe3, Ue3)
    compare_step_ele(dm1, Xe1, Ue1, dm3, Xe3, Ue3)
    #compare_step_ailerons(dm1, Xe1, Ue1, dm3, Xe3, Ue3)
    #compare_step_rudder(dm1, Xe1, Ue1, dm3, Xe3, Ue3)
    #compare_step_flaps(dm1, Xe1, Ue1, dm3, Xe3, Ue3)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    np.set_printoptions(linewidth=500)
    main()
