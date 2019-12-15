import math, numpy as np

import pat3.algebra as p3_alg

import pdb

def R_aero_to_body(alpha, beta):
    """
    computes the aero to body rotation matrix
    """
    ca, sa = np.cos(alpha), np.sin(alpha) 
    cb, sb = np.cos(beta), np.sin(beta)
    return np.array([[ca*cb, -ca*sb, -sa],
                     [sb   ,  cb   ,  0.],
                     [sa*cb, -sa*sb,  ca]])

def vel_world_to_aero_eul(ivel_world, eul, wind_ned):
    '''
    '''
    world_to_body_R = p3_alg.rmat_of_euler(eul)
    return vel_world_to_aero_R(ivel_world, world_to_body_R, wind_ned)

def vel_world_to_aero_quat(ivel_world, quat, wind_ned):
    '''
    '''
    world_to_body_R = p3_alg.rmat_of_quat(quat)
    return vel_world_to_aero_R(ivel_world, world_to_body_R, wind_ned)

def vel_world_to_aero_R(ivel_world, world_to_body_R, wind_ned):
    avel_world = np.asarray(ivel_world) - wind_ned
    va = np.linalg.norm(avel_world)
    avel_body = u, v, w = np.dot(world_to_body_R, avel_world)
    alpha = np.arctan2(w, u)
    beta = np.arcsin(v/va)
    return va, alpha, beta
    


def vel_aero_to_world(avel_aero, eul, wind_ned):
    '''
    '''
    va, alpha, beta = avel_aero
    _avel_aero = np.array([va, 0, 0])
    avel_body = np.dot(R_aero_to_body(alpha, beta), _avel_aero)
    body_to_world_R = p3_alg.rmat_of_euler(eul).T
    avel_world = np.dot(body_to_world_R, avel_body)
    ivel_world = wind_ned + avel_world
    return ivel_world


class SixDOFEuclidianEuler():
    sv_x, sv_y, sv_z, sv_xd, sv_yd, sv_zd, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)

    sv_slice_pos   = slice(sv_x,   sv_z+1)    # position ned (north east down)
    sv_slice_vel   = slice(sv_xd,   sv_zd+1)  # inertial velocity in world frame (ned)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)  # euler angles from world to body ( FIXME check direction )
    sv_slice_rvel  = slice(sv_p,   sv_r+1)    # rotational velocity in body frame

    @classmethod
    def to_six_dof_aero_euler(cls, Xee, atm=None, t=0):
        Xae = np.zeros(SixDOFAeroEuler.sv_size)
        Xae[SixDOFAeroEuler.sv_slice_pos]   = Xee[cls.sv_slice_pos]
        Xae[SixDOFAeroEuler.sv_slice_eul]   = Xee[cls.sv_slice_eul]
        Xae[SixDOFAeroEuler.sv_slice_rvel]  = Xee[cls.sv_slice_rvel]
        wind_ned = (atm.get_wind(Xee[cls.sv_slice_pos], t) if atm is not None else [0, 0, 0])
        Xae[SixDOFAeroEuler.sv_slice_vaero] = vel_world_to_aero(Xee[cls.sv_slice_vel], Xee[cls.sv_slice_eul], wind_ned)
        return Xae

    @classmethod
    def from_six_dof_aero_euler(cls, Xae, atm=None, t=0):
        return SixDOFAeroEuler.to_six_dof_euclidian_euler(Xae)

    @classmethod
    def cont_dyn(cls, X, t, U, P, add_weight=False):
        Fb, Mb = U[:3], U[3:]
        Xd = np.zeros(cls.sv_size)
        p_w, v_w, e_w2b, om_b = X[cls.sv_slice_pos], X[cls.sv_slice_vel], X[cls.sv_slice_eul], X[cls.sv_slice_rvel]
        # Translational kinematics
        Xd[cls.sv_slice_pos] = v_w
        # Newton for forces
        R_b2w =  p3_alg.rmat_of_euler(e_w2b).T
        Xd[cls.sv_slice_vel] = 1./P.m*(np.dot(R_b2w, Fb))
        # do we add weight or not???
        if add_weight: Xd[cls.sv_slice_vel] += 1./P.m*(np.array([0, 0, P.m*P.g]))
        # Rotational kinematics
        Xd[cls.sv_slice_eul] = p3_alg.euler_derivatives(e_w2b, om_b)
        # Newton for moments
        Xd[cls.sv_slice_rvel] = np.dot(P.invI, Mb - np.cross(om_b, np.dot(P.I, om_b)))
        return Xd

        

class SixDOFAeroEuler():
    sv_x, sv_y, sv_z, sv_va, sv_alpha, sv_beta, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)
    sv_slice_pos   = slice(sv_x,   sv_z+1)
    sv_slice_vaero = slice(sv_va,  sv_beta+1)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)
    sv_slice_rvel  = slice(sv_p,   sv_r+1)

    @classmethod
    def to_six_dof_euclidian_euler(cls, X, atm=None, t=0):
        Xee = np.zeros(SixDOFEuclidianEuler.sv_size)
        Xee[SixDOFEuclidianEuler.sv_slice_pos] = X[cls.sv_slice_pos]
        Xee[SixDOFEuclidianEuler.sv_slice_eul] = X[cls.sv_slice_eul]
        Xee[SixDOFEuclidianEuler.sv_slice_rvel] = X[cls.sv_slice_rvel]
        wind_ned = (atm.get_wind(X[cls.sv_slice_pos], t) if atm is not None else [0, 0, 0])
        Xee[SixDOFEuclidianEuler.sv_slice_vel] = vel_aero_to_world(X[cls.sv_slice_vaero], X[cls.sv_slice_eul], wind_ned)
        return Xee

    @classmethod
    def to_six_dof_euclidian_quat(cls, Xae, atm=None, t=0):
        Xeq = np.zeros(SixDOFEuclidianQuat.sv_size)
        Xeq[SixDOFEuclidianQuat.sv_slice_pos] = Xae[cls.sv_slice_pos]
        wind_ned = (atm.get_wind(Xae[cls.sv_slice_pos], t) if atm is not None else [0, 0, 0])
        Xeq[SixDOFEuclidianEuler.sv_slice_vel] = vel_aero_to_world(Xae[cls.sv_slice_vaero], Xae[cls.sv_slice_eul], wind_ned)
        Xeq[SixDOFEuclidianQuat.sv_slice_quat] = p3_alg.quat_of_euler(Xae[cls.sv_slice_eul])
        Xeq[SixDOFEuclidianEuler.sv_slice_rvel] = Xae[cls.sv_slice_rvel]
        return Xeq

    @classmethod
    def from_six_dof_aero_euler(cls, X): return X
    
    @classmethod
    def state_str(cls, X):
        return """pos: {:-.2f}, {:-.2f}, {:-.2f} m
        vel: {:-.2f} m/s, alpha {:-.2f}, beta {:-.2f} deg
        att:    {:-.2f}, {:-.2f}, {:-.2f} deg
        """.format(X[cls.sv_x], X[cls.sv_y], X[cls.sv_z],
                   X[cls.sv_va], np.rad2deg(X[cls.sv_alpha]), np.rad2deg(X[cls.sv_beta]),
                   np.rad2deg(X[cls.sv_phi]), np.rad2deg(X[cls.sv_theta]), np.rad2deg(X[cls.sv_psi]))

    @classmethod
    def cont_dyn(cls, X, t, U, P, atm=None):
        Fb, Mb = U[:3], U[3:]
        X_pos = X[cls.sv_slice_pos]                      # ned pos
        X_avel = va, alpha, beta = X[cls.sv_slice_vaero] # airvel, alpha, beta
        X_euler = X[cls.sv_slice_eul]                    # euler angles
        X_rvel_body = X[cls.sv_slice_rvel]               # body rotational velocities

        earth_to_body_R = p3_alg.rmat_of_euler(X_euler)
        p_w, v_aero, e_w2b, om_b = X[cls.sv_slice_pos], X[cls.sv_slice_vaero], X[cls.sv_slice_eul], X[cls.sv_slice_rvel]
        wind_ned = atm.get_wind(X_pos, t) if atm is not None else [0, 0, 0]
        waccel_body = [0, 0, 0]  # np.dot(earth_to_body, waccel_earth)
        ivel_world = vel_aero_to_world(X_avel, X_euler, wind_ned)
        ivel_body = np.dot(earth_to_body_R, ivel_world)
        avel_aero = [X[cls.sv_va], 0., 0.]
        _R_aero_to_body = R_aero_to_body(alpha, beta)
        avel_body = np.dot(_R_aero_to_body, avel_aero)
        
        Xd = np.zeros(cls.sv_size)
        # Translational kinematics
        Xd[cls.sv_slice_pos] = vel_aero_to_world(X_avel, X_euler, wind_ned)
        # Translational dynamics
        iaccel_body = 1./P.m*Fb - np.cross(X_rvel_body, ivel_body)
        aaccel_body = iaccel_body - waccel_body
        Xd[cls.sv_va] = np.inner(avel_body, aaccel_body)/X[cls.sv_va]
        (avx, avy, avz), (aax, aay, aaz) = avel_body, aaccel_body
        Xd[cls.sv_alpha] = (avx*aaz - avz*aax)/(avx**2+avz**2)
        Xd[cls.sv_beta] = (X[cls.sv_va]*aay - avy*Xd[cls.sv_va]) / X[cls.sv_va] / math.sqrt(avx**2+aaz**2)
        # Euler angles kinematics
        Xd[cls.sv_slice_eul] = p3_alg.euler_derivatives(X_euler, X_rvel_body)
        # Rotational dynamics
        Xd[cls.sv_slice_rvel] = np.dot(P.invI, Mb - np.cross(X_rvel_body, np.dot(P.I, X_rvel_body)))

        return Xd
        #return np.zeros(SixDOFAeroEuler.sv_size)

class SixDOFBodyEuler():
    
    sv_x, sv_y, sv_z, sv_u, sv_v, sv_w, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)

    sv_slice_pos   = slice(sv_x,   sv_z+1)    # position ned (north east down)
    sv_slice_vel   = slice(sv_u,   sv_w+1)    # inertial velocity in body frame (frd)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)  # euler angles
    sv_slice_rvel  = slice(sv_p,   sv_r+1)    # rotational velocity in body frame
    

class SixDOFEuclidianQuat():
    
    sv_x, sv_y, sv_z, sv_xd, sv_yd, sv_zd, sv_qi, sv_qx, sv_qy, sv_qz, sv_p, sv_q, sv_r, sv_size = range(0,14)
    
    sv_slice_pos   = slice(sv_x,  sv_z+1)   # position ned (north east down)
    sv_slice_vel   = slice(sv_xd, sv_zd+1)  # inertial velocity in world frame (ned)
    sv_slice_quat  = slice(sv_qi, sv_qz+1)  # quaternion from world to body ( FIXME check direction )
    sv_slice_rvel  = slice(sv_p,  sv_r+1)   # rotational velocity in body frame

    @classmethod
    def from_six_dof_aero_euler(cls, Xae, atm=None, t=0):
        return SixDOFAeroEuler.to_six_dof_euclidian_quat(Xae, atm, t)

    @classmethod
    def to_six_dof_euclidian_euler(cls, Xeq, atm=None, t=0):
        Xee = np.zeros(SixDOFEuclidianEuler.sv_size)
        Xee[SixDOFEuclidianEuler.sv_slice_pos] = Xeq[cls.sv_slice_pos]
        Xee[SixDOFEuclidianEuler.sv_slice_vel] = Xeq[cls.sv_slice_vel]
        Xee[SixDOFEuclidianEuler.sv_slice_eul] = p3_alg.euler_of_quat(Xeq[cls.sv_slice_quat])
        Xee[SixDOFEuclidianEuler.sv_slice_rvel] = Xeq[cls.sv_slice_rvel]
        return Xee
    
    @classmethod
    def to_six_dof_aero_euler(cls, Xeq, atm=None, t=0):
        Xae = np.zeros(SixDOFAeroEuler.sv_size)
        Xae[SixDOFAeroEuler.sv_slice_pos]   = Xeq[cls.sv_slice_pos]
        wind_ned = (atm.get_wind(Xeq[cls.sv_slice_pos], t) if atm is not None else [0, 0, 0])
        Xae[SixDOFAeroEuler.sv_slice_vaero] = vel_world_to_aero(Xeq[cls.sv_slice_vel],  Xae[SixDOFAeroEuler.sv_slice_eul], wind_ned)
        Xae[SixDOFAeroEuler.sv_slice_eul]   = p3_alg.euler_of_quat(Xeq[cls.sv_slice_quat])
        Xae[SixDOFAeroEuler.sv_slice_rvel]  = Xeq[cls.sv_slice_rvel]
        return Xae
    
    @classmethod
    def cont_dyn(cls, X, t, U, P, add_weight=False):
        Fb, Mb = U[:3], U[3:]
        Xd = np.zeros(cls.sv_size)
        p_w, v_w, q_w2b, om_b = X[cls.sv_slice_pos], X[cls.sv_slice_vel], X[cls.sv_slice_quat], X[cls.sv_slice_rvel]
        # Translational kinematics
        Xd[cls.sv_slice_pos] = v_w
        # Newton for forces
        R_b2w =  p3_alg.rmat_of_quat(q_w2b).T
        Xd[cls.sv_slice_vel] = 1./P.m*(np.dot(R_b2w, Fb))
        # do we add weight or not???
        if add_weight: Xd[cls.sv_slice_vel] += 1./P.m*(np.array([0, 0, P.m*P.g]))
        # Rotational kinematics
        Xd[cls.sv_slice_quat] = p3_alg.quat_derivative(q_w2b, om_b)
        # Newton for moments
        Xd[cls.sv_slice_rvel] = np.dot(P.invI, Mb - np.cross(om_b, np.dot(P.I, om_b)))
        return Xd
