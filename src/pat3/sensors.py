#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#
# Copyright 2013-2021 Antoine Drouin (poinix@gmail.com)
#
# This file is part of PAT.
#
#    PAT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PAT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PAT.  If not, see <http://www.gnu.org/licenses/>.
#

import pdb

"""
  sensors related stuff
"""

import math, numpy as np, matplotlib.pyplot as plt
from numpy.random import default_rng
import scipy.interpolate

import pat3.frames as p3_fr

class Sensors:
    def __init__(self, std_va=0.1):
        self.va_std = std_va
        self.rng = default_rng()

    def get_measurements(self, Xae, Xee):
        Xae2 = np.array(Xae)
        Xae2[p3_fr.SixDOFAeroEuler.sv_va] +=  self.va_std * self.rng.standard_normal(1)
        return Xae2, Xee




