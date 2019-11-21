import numpy as np

import pat3.algebra as p3_alg

import pdb

def R_aero_to_body(alpha, beta):
    """
    computes the aero to body rotation matix
    """
    ca, sa = np.cos(alpha), np.sin(alpha) 
    cb, sb = np.cos(beta), np.sin(beta)
    return np.array([[ca*cb, -ca*sb, -sa],
                     [sb   ,  cb   ,  0.],
                     [sa*cb, -sa*sb,  ca]])


class SixDOFAeroEuler():
    sv_x, zv_y, sv_z, sv_va, sv_alpha, sv_beta, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)
    sv_slice_pos   = slice(sv_x,   sv_z+1)
    sv_slice_vaero = slice(sv_va,   sv_beta+1)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)
    sv_slice_rvel  = slice(sv_p,   sv_r+1)

    @classmethod
    def to_six_dof_euclidian_euler(cls, X, atm=None):
        Xee = np.zeros(SixDOFEuclidianEuler.sv_size)
        Xee[SixDOFEuclidianEuler.sv_slice_pos] = X[cls.sv_slice_pos]
        Xee[SixDOFEuclidianEuler.sv_slice_eul] = X[cls.sv_slice_eul]
        Xee[SixDOFEuclidianEuler.sv_slice_rvel] = X[cls.sv_slice_rvel]
        alpha, beta = X[cls.sv_alpha], X[cls.sv_alpha]
        aero_to_body = R_aero_to_body(alpha, beta)
        body_to_earth = p3_alg.rmat_of_euler(X[cls.sv_slice_eul]).T
        avel_aero = [X[cls.sv_va], 0., 0.]
        avel_body = np.dot(aero_to_body, avel_aero)
        avel_earth = np.dot(body_to_earth, avel_body)
        wvel_earth = (atm.get_wind(X[cls.sv_slice_pos], t=0) if atm is not None else [0, 0, 0])
        Xee[SixDOFEuclidianEuler.sv_slice_v] = avel_earth + wvel_earth
        return Xee
    
class SixDOFEuclidianEuler():
    sv_x, sv_y, sv_z, sv_xd, sv_yd, sv_zd, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)

    sv_slice_pos   = slice(sv_x,   sv_z+1)
    sv_slice_v     = slice(sv_xd,   sv_zd+1)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)
    sv_slice_rvel  = slice(sv_p,   sv_r+1)



