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
# 
def vel_world_to_aero(pos_ned, ivel_world, eul, atm=None):
    wvel_world = (atm.get_wind(pos_ned, t=0) if atm is not None else [0, 0, 0])
    avel_world = np.asarray(ivel_world) - wvel_world
    world_to_body_R = p3_alg.rmat_of_euler(eul)
    va = np.linalg.norm(avel_world)
    avel_body = u, v, w = np.dot(world_to_body_R, avel_world)
    alpha = np.arctan2(w, u)
    beta = np.arctan2(v, va)
    return va, alpha, beta

#
def vel_aero_to_world(pos_ned, avel_aero, eul, atm=None):
    avel_body = u, v, w = np.dot(R_aero_to_body(avel_aero[1], avel_aero[2]), avel_aero)
    body_to_world_R = p3_alg.rmat_of_euler(eul).T
    avel_world = np.dot(body_to_world_R, avel_body)
    wvel_world = (atm.get_wind(pos_ned, t=0) if atm is not None else [0, 0, 0])
    ivel_world = wvel_world + avel_world
    return ivel_world


class SixDOFAeroEuler():
    sv_x, sv_y, sv_z, sv_va, sv_alpha, sv_beta, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)
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
        aero_to_body = R_aero_to_body(X[cls.sv_alpha], X[cls.sv_beta])
        body_to_earth = p3_alg.rmat_of_euler(X[cls.sv_slice_eul]).T
        avel_aero = [X[cls.sv_va], 0., 0.]
        avel_body = np.dot(aero_to_body, avel_aero)
        avel_earth = np.dot(body_to_earth, avel_body)
        wvel_earth = (atm.get_wind(X[cls.sv_slice_pos], t=0) if atm is not None else [0, 0, 0])
        Xee[SixDOFEuclidianEuler.sv_slice_vel] = avel_earth + wvel_earth
        return Xee

    @classmethod
    def state_str(cls, X):
        return """pos: {:-.2f}, {:-.2f}, {:-.2f} m
        vel: {:-.2f} m/s, alpha {:-.2f}, beta {:-.2f} deg
        att:    {:-.2f}, {:-.2f}, {:-.2f} deg
        """.format(X[cls.sv_x], X[cls.sv_y], X[cls.sv_z],
                   X[cls.sv_va], np.rad2deg(X[cls.sv_alpha]), np.rad2deg(X[cls.sv_beta]),
                   np.rad2deg(X[cls.sv_phi]), np.rad2deg(X[cls.sv_theta]), np.rad2deg(X[cls.sv_psi]))

    
class SixDOFEuclidianEuler():
    sv_x, sv_y, sv_z, sv_xd, sv_yd, sv_zd, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)

    sv_slice_pos   = slice(sv_x,   sv_z+1)    # position ned (north east down)
    sv_slice_vel   = slice(sv_xd,   sv_zd+1)  # inertial velocity in world frame (ned)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)  # euler angles
    sv_slice_rvel  = slice(sv_p,   sv_r+1)    # rotational velocity in body frame

    @classmethod
    def to_six_dof_aero_euler(cls, Xee, atm=None):
        Xae = np.zeros(SixDOFAeroEuler.sv_size)
        Xae[SixDOFAeroEuler.sv_slice_pos]   = Xee[cls.sv_slice_pos]
        Xae[SixDOFAeroEuler.sv_slice_eul]   = Xee[cls.sv_slice_eul]
        Xae[SixDOFAeroEuler.sv_slice_rvel]  = Xee[cls.sv_slice_rvel]
        Xae[SixDOFAeroEuler.sv_slice_vaero] = vel_world_to_aero(Xee[cls.sv_slice_pos], Xee[cls.sv_slice_vel], Xee[cls.sv_slice_eul], atm)
        return Xae
