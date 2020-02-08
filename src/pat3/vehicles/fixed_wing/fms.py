import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt

import pat3.vehicles.fixed_wing.legacy_6dof as p1_fw_dyn
import pat3.vehicles.fixed_wing.piloting as p3_pil

import pat3.utils as p3_u
import pat3.algebra as p3_alg
import pat3.frames as p3_fr
import pat3.trajectory_3D as p3_traj3d
import pat3.plot_utils as p3_pu

import pat3.vehicles.fixed_wing.guidance as p3_guid
import pat3.vehicles.fixed_wing.guidance_ardusoaring as p3_guidardu
import pat3.vehicles.fixed_wing.guidance_soaring as p3_guidsoar

class FMS:
    mod_auto1, mod_circle, mod_centering, mod_ardu, mod_searching, mod_nb = range(6)
    def __init__(self, dm, trim_args, dt=0.01):
        self.mode = FMS.mod_circle
        self.Xe, self.Ue = dm.trim(trim_args, report=True)
        self.guidances = [ p3_guid.GuidanceAuto1(self.Xe, self.Ue, dm, dt),
                           p3_guid.GuidanceCircle(dm, trim_args, dt=0.01),
                           p3_guidsoar.GuidanceSoaring(dm, trim_args),
                           p3_guidardu.GuidanceArduSoaring(dm, trim_args),
                           p3_guid.GuidanceSearching(dm, trim_args)]
        self.spv_size = 2
        self.last_X, self.last_t  = self.Xe, 0 # FIXME, move that elsewhere?
    
    def set_mode(self, m):
        if self.mode == FMS.mod_circle and m == FMS.mod_centering:
            #self.thermaling_traj.c = self.traj.c
            #self.thermaling_traj.r = self.traj.r
            self.guidances[FMS.mod_centering].set_circle_center(self.guidances[FMS.mod_circle].traj.c)

        self.mode = m
        self.guidances[self.mode].enter(self.last_X, self.last_t) # t, X ?

    def guidance(self): return self.guidances[self.mode]
        
    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        self.last_X, self.last_t = X, t
        return self.guidances[self.mode].get(t, X, Xee, Yc, debug, traj)

    def carrot(self): return self.guidance().carrot
    def phi_sp(self): return self.guidance().phi_sp
    def theta_sp(self): return self.guidance().theta_sp
    def get_va_sp(self): return self.guidance().v_sp
    def set_va_sp(self, _v): self.guidance().v_sp = _v
