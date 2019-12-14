import sys, os, math, numpy as np
import logging
import pdb

import matplotlib.pyplot as plt

import pat3.vehicles.fixed_wing.guidance as p3_guid
import pat3.trajectory_3D as p3_traj3d


   
class GuidanceArduSoaring(p3_guid.GuidancePurePursuit):
    def __init__(self, dm, traj, trim_args = {'h':0, 'va':12, 'gamma':0}, dt=0.01):
        p3_guid.GuidancePurePursuit.__init__(self, dm, p3_traj3d.CircleRefTraj(c=[0., 15., 0.], r=15.), trim_args, dt)


    def get(self, t, X, Xee=None, Yc=None, debug=False, traj=None):
        return p3_guid.GuidancePurePursuit.get(self, t, X, Xee, Yc, debug)
