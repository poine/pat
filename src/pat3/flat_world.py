#-*- coding: utf-8 -*-
#
# Copyright 2013-2015 Antoine Drouin (poinix@gmail.com)
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

"""
 This modules handles frame convertions
 
"""

import math, numpy as np



""" Rotation Matrices """
def dcm_eci_to_ecef(utc):
    ''' 
    Returns the rotation matrix from ECI (Earth Centered Inertial) to ECEF (Earth Centered Earth Fixed)
      Uses IAU-2000/2006 (International Astronomical Union (IAU)-2000/2005 reference system)

    utc: UTC (Universal Coordinated Time) as a datetime object

    see: http://www.iers.org
    http://www.mathworks.com/matlabcentral/fileexchange/28233-convert-eci-to-ecef-coordinates
    '''
    return np.identity(3)



def earth_to_body((phi, theta, psi)):
    """ computes the earth to body rotation matix """
    cphi, sphi     = math.cos(phi),   math.sin(phi)
    ctheta, stheta = math.cos(theta), math.sin(theta)
    cpsi, spsi     = math.cos(psi),   math.sin(psi)
    return np.array([[ctheta*cpsi                 , ctheta*spsi                 , -stheta],
                     [sphi*stheta*cpsi - cphi*spsi, sphi*stheta*spsi + cphi*cpsi, sphi*ctheta],
                     [cphi*stheta*cpsi + sphi*spsi, cphi*stheta*spsi - sphi*cpsi, cphi*ctheta]])


def body_to_earth((phi, theta, psi)): return earth_to_body((phi, theta, psi)).T


def aero_to_body(alpha, beta):
    """ computes the aero to body rotation matix """
    ca = math.cos(alpha); sa = math.sin(alpha) 
    cb = math.cos(beta);  sb = math.sin(beta)
    return np.array([[ca*cb, -ca*sb, -sa],
                     [sb   ,  cb   ,  0.],
                     [sa*cb, -sa*sb,  ca]])


def body_to_aero(alpha, beta): return aero_to_body(alpha, beta).T

def aero_to_earth(eul, (alpha, beta)):
    return np.dot(body_to_earth(eul), aero_to_body(alpha, beta))



""" Rotation Matrices kinematics """
def earth_to_body_dot(eul, eul_dot):
    ph, th, ps = eul
    phd, thd, psd = eul_dot
    cph, sph = math.cos(ph), math.sin(ph)
    cth, sth = math.cos(th), math.sin(th)
    cps, sps = math.cos(ps), math.sin(ps)
    a1 = -thd*sth*cps-psd*cth*sps
    a2 = -thd*sth*sps+psd*cth*cps
    a3 = -thd*cth
    a4 = phd*(cph*sth*cps+sph*sps)+thd*(sph*cth*cps)+psd*(-sph*sth*sps-cph*cps)
    a5 = phd*(cph*sth*sps-sph*cps)+thd*(sph*cth*sps)+psd*( sph*sth*cps-cph*sps)
    a6 = phd*cph*cth - thd*sph*sth
    a7 = phd*(-sph*sth*cps+cph*sps)+thd*(cph*cth*cps)+psd*(-cph*sth*sps+sph*cps)
    a8 = phd*(-sph*sth*sps-cph*cps)+thd*(cph*cth*sps)+psd*( cph*sth*cps+sph*sps)
    a9 = -phd*sph*cth - thd*cph*sth
    return np.array([[a1, a2, a3],[a4, a5, a6],[a7, a8, a9]])


def body_to_earth_dot(eul, eul_dot):
    return earth_to_body_dot(eul, eul_dot).T


def aero_to_body_dot(alpha, beta, alpha_dot, beta_dot):
    ca, sa, cb, sb = math.cos(alpha), math.sin(alpha), math.cos(beta), math.sin(beta)
    cacb, casb, sacb, sasb = ca*cb, ca*sb, sa*cb, sa*sb
    return np.array([[-alpha_dot*sacb-beta_dot*casb, alpha_dot*sasb-beta_dot*cacb,  -alpha_dot*ca],
                     [beta_dot*cb,                   -beta_dot*sb,                   0.],
                     [ alpha_dot*cacb-beta_dot*sasb, -alpha_dot*casb-beta_dot*sacb, -alpha_dot*sa]])
    


import pdb

def aero_to_earth_dot(vi_e, w_e, eul, aci_e, eul_d):
    b2e = body_to_earth(eul)
    b2e_dot = body_to_earth_dot(eul, eul_d)
    va, alpha, beta, va_dot, alpha_dot, beta_dot = get_va_alpha_beta_dot(vi_e, w_e, eul, aci_e, eul_d)
    a2b = aero_to_body(alpha, beta)
    a2b_dot = aero_to_body_dot(alpha, beta, alpha_dot, beta_dot)
    return np.dot(b2e_dot, a2b) + np.dot(b2e, a2b_dot)



#

def get_va_alpha_beta(vi_e, w_e, eul):
    '''
    vi_e:  inertial velocity in earth frame
    w_e:   wind in earth frame
    eul:   euler angles
    '''
    va_e = vi_e - w_e
    va = np.linalg.norm(va_e)
    e2b = earth_to_body(eul)
    va_b = np.dot(e2b, va_e)
    alpha, beta = math.atan(va_b[2]/va_b[0]), math.asin(va_b[1]/va)
    return va, alpha, beta 


def get_va_alpha_beta_dot(vi_e, w_e, eul, aci_e, eul_d):
    '''
    vi_e:  inertial velocity in earth frame
    w_e:   wind in earth frame
    eul:   euler angles
    aci_e: inertial acceleration in earth frame
    eul_d: time derivative of euler angles
    '''
    va_e = vi_e - w_e
    va = np.linalg.norm(va_e)
    e2b = earth_to_body(eul)
    va_b = np.dot(e2b, va_e)
    f = va_b[2]/va_b[0]
    alpha, beta = math.atan(f), math.asin(va_b[1]/va) 
    va_dot = np.inner(aci_e, va_e)/va
    
    e2b_dot =  earth_to_body_dot(eul, eul_d)
    va_b_dot = np.dot(e2b_dot, va_e) + np.dot(e2b, aci_e)
    
    f_dot = (va_b_dot[2]*va_b[0]-va_b[2]*va_b_dot[0])/(va_b[0]**2)
    alpha_dot = f_dot/(1+f**2)

    f = va_b[1]/va
    f_dot = (va_b_dot[1]*va-va_b[1]*va_dot)/(va**2)
    beta_dot = f_dot/math.sqrt(1-f**2)

    return va, alpha, beta, va_dot, alpha_dot, beta_dot
