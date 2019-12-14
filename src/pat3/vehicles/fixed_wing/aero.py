#-*- coding: utf-8 -*-
#
# Copyright 2019 Antoine Drouin (poinix@gmail.com)
#
# This file is part of PAT.
#
#    PAT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PAT is distributed in the hpuope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with PAT.  If not, see <http://www.gnu.org/licenses/>.
#
import os
import math, numpy as np

import pdb
import pat3.frames as p3_fr
import pat3.atmosphere as p3_atm

'''
 Aerodynamics for fixed wing aircraft
'''


#
# Forces
#

def thrust_of_throttle(throttle, fmax, rho, rhor, nrho, v, vr, nv):
    return throttle*fmax*math.pow((rho/rhor),nrho)*math.pow((v/vr),nv)

def get_f_eng_body(h, va, U, P):
    """
    return propulsion forces expressed in body frame
    """
    rho = p3_atm.get_rho(h)
    thrusts = [thrust_of_throttle(U[i], P.fmaxs[i], rho, P.rhois[i], P.nrhos[i], va, P.Vis[i], P.nVs[i]) for i in range(P.eng_nb)]
    f_engines_body = np.array([np.dot(P.eng_to_body[i], np.array([thrusts[i], 0., 0.])) for i in range(P.eng_nb)])
    return f_engines_body 

def get_f_aero_coef(alpha, beta, rvel, Usfc, P):
    """
    return aero coefficients for forces
    """
    d_alpha = alpha - P.alpha0
    nrvel = rvel*np.array([P.Bref, P.Cref, P.Bref])/2/P.Vref # FIXME va??
    CL = P.CL0 + P.CL_alpha*d_alpha + P.CL_beta*beta + np.dot(P.CL_omega,nrvel) + np.dot(P.CL_sfc,Usfc)
    CD = P.CD0 + P.CD_k1*CL + P.CD_k2*(CL**2) + np.dot(P.CD_sfc,Usfc)
    CY = P.CY_alpha*d_alpha + P.CY_beta*beta + np.dot(P.CY_omega,nrvel) + np.dot(P.CY_sfc,Usfc)
    return [CL, CY, CD]

def get_f_aero_body(va, alpha, beta, rvel, Usfc, P, Pdyn):
    """
    return aerodynamic forces in body frame
    """
    CL, CY, CD = get_f_aero_coef(alpha, beta, rvel, Usfc, P)
    F_aero_body = Pdyn*P.Sref*np.dot(p3_fr.R_aero_to_body(alpha, beta), [-CD, CY, -CL])
    return F_aero_body


#
# Moments
#

def get_m_eng_body(f_eng_body, P):
    """
    return propulsion moments expressed in body frame
    """
    m = np.zeros(3)
    for i in range(0, P.eng_nb):
        m += np.cross(P.eng_pos[i], f_eng_body[i])
    return m

def get_m_aero_coef(alpha, beta, rvel, Usfc, P):
    d_alpha = alpha - P.alpha0
    nrvel =  rvel*np.array([P.Bref, P.Cref, P.Bref])/2/P.Vref
    Cl =         P.Cl_alpha*d_alpha + P.Cl_beta*beta +\
         np.dot(P.Cl_omega,nrvel) + np.dot(P.Cl_sfc,Usfc)
    Cm = P.Cm0 + P.Cm_alpha*d_alpha + P.Cm_beta*beta +\
         np.dot(P.Cm_omega,nrvel) + np.dot(P.Cm_sfc,Usfc)
    Cn =         P.Cn_alpha*d_alpha + P.Cn_beta*beta +\
         np.dot(P.Cn_omega,nrvel) + np.dot(P.Cn_sfc,Usfc)
    return Cl, Cm, Cn

def get_m_aero_body(va, alpha, beta, rvel, Usfc, P, Pdyn):
    Cl, Cm, Cn = get_m_aero_coef(alpha, beta, rvel, Usfc, P)
    return Pdyn*P.Sref*np.array([Cl*P.Bref, Cm*P.Cref, Cn*P.Bref])
