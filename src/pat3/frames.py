import math, numpy as np
import matplotlib.pyplot as plt

import pat3.algebra as p3_alg
import pat3.plot_utils as p3_pu

import pdb

def ned_to_enu(X):
    return np.vstack([X[:,1], X[:,0], -X[:,2]]).T


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
    '''
    '''
    avel_world = np.asarray(ivel_world) - wind_ned
    va = np.linalg.norm(avel_world)
    avel_body = u, v, w = np.dot(world_to_body_R, avel_world)
    alpha = np.arctan2(w, u)
    beta = np.arcsin(v/va)
    return va, alpha, beta
    
def vel_aero_to_world_euler(avel_aero, eul, wind_ned):
    '''
    '''
    body_to_world_R = p3_alg.rmat_of_euler(eul).T
    return vel_aero_to_world_R(avel_aero, body_to_world_R, wind_ned)

def vel_aero_to_world_R(avel_aero, body_to_world_R, wind_ned):
    '''
    '''
    va, alpha, beta = avel_aero
    _avel_aero = np.array([va, 0, 0])
    avel_body = np.dot(R_aero_to_body(alpha, beta), _avel_aero)
    avel_world = np.dot(body_to_world_R, avel_body)
    ivel_world = wind_ned + avel_world
    return ivel_world

def vel_aero_to_body_R(avel_aero, world_to_body_R, wind_ned):
    '''
    '''
    va, alpha, beta = avel_aero
    _avel_aero = np.array([va, 0, 0])
    avel_body = np.dot(R_aero_to_body(alpha, beta), _avel_aero)
    wind_body = np.dot(world_to_body_R, wind_ned)
    ivel_body = wind_body + avel_body
    return ivel_body


def vel_world_to_body_R(ivel_world, world_to_body_R):
    return np.dot(world_to_body_R, ivel_world)

def vel_world_to_body_eul(ivel_world, eul):
    world_to_body_R = p3_alg.rmat_of_euler(eul)
    return vel_world_to_body_R(ivel_world, world_to_body_R)    

def vel_body_to_world_R(ivel_body, body_to_world_R):
    return np.dot(body_to_world_R, ivel_body)

def vel_body_to_world_euler(ivel_body, eul):
    body_to_world_R = p3_alg.rmat_of_euler(eul).T
    return vel_body_to_world_R(ivel_body, body_to_world_R)

def vel_body_to_world_quat(ivel_body, quat):
    body_to_world_R = p3_alg.rmat_of_quat(quat).T
    return vel_body_to_world_R(ivel_body, body_to_world_R)


def plot_input(time, U, figure):
    Fb, Mb = U[:,:3], U[:,3:]
    ax = figure.add_subplot(5, 3, 13)
    plt.plot(time, Fb[:,0], label='Fx')
    plt.plot(time, Fb[:,1], label='Fy')
    plt.plot(time, Fb[:,2], label='Fz')
    p3_pu.decorate(ax, title="$F_b$", ylab="N", min_yspan=0.001, legend=True)
    ax = figure.add_subplot(5, 3, 14)
    plt.plot(time, Mb[:,0], label='Mx')
    plt.plot(time, Mb[:,1], label='My')
    plt.plot(time, Mb[:,2], label='Mz')
    p3_pu.decorate(ax, title="$M_b$", ylab="Nm", min_yspan=0.001, legend=True)
    

#
#
#   Euclidian/Euler
#
#
class SixDOFEuclidianEuler():
    sv_x, sv_y, sv_z, sv_xd, sv_yd, sv_zd, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)

    sv_slice_pos   = slice(sv_x,   sv_z+1)    # position ned (north east down)
    sv_slice_vel   = slice(sv_xd,   sv_zd+1)  # inertial velocity in world frame (ned)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)  # euler angles from world to body ( FIXME check direction )
    sv_slice_rvel  = slice(sv_p,   sv_r+1)    # rotational velocity in body frame

    @classmethod
    def get_name(cls): return 'euclidian_euler'
    @classmethod
    def get_default_state(cls): return np.zeros(cls.sv_size)
    
    @classmethod
    def to_six_dof_euclidian_quat(cls, Xee, atm=None, t=0):
        Xeq = np.zeros(SixDOFEuclidianQuat.sv_size)
        Xeq[SixDOFEuclidianQuat.sv_slice_pos] = Xee[cls.sv_slice_pos]
        Xeq[SixDOFEuclidianQuat.sv_slice_vel] = Xee[cls.sv_slice_vel]
        Xeq[SixDOFEuclidianQuat.sv_slice_quat] = p3_alg.quat_of_euler(Xee[cls.sv_slice_eul])
        Xeq[SixDOFEuclidianQuat.sv_slice_rvel] = Xee[cls.sv_slice_rvel]
        return Xeq

    @classmethod
    def to_six_dof_body_euler(cls, Xee, atm=None, t=0):
        Xbe = np.zeros(SixDOFBodyEuler.sv_size)
        Xbe[SixDOFBodyEuler.sv_slice_pos] = Xee[cls.sv_slice_pos]
        Xbe[SixDOFBodyEuler.sv_slice_vel] = vel_world_to_body_eul(Xee[cls.sv_slice_vel], Xee[cls.sv_slice_eul])
        Xbe[SixDOFBodyEuler.sv_slice_eul] = Xee[cls.sv_slice_eul]
        Xbe[SixDOFBodyEuler.sv_slice_rvel] = Xee[cls.sv_slice_rvel]
        return Xbe

    @classmethod
    def to_six_dof_body_quat(cls, Xee, atm=None, t=0):
        Xbq = np.zeros(SixDOFBodyQuat.sv_size)
        Xbq[SixDOFBodyQuat.sv_slice_pos] = Xee[cls.sv_slice_pos]
        Xbq[SixDOFBodyQuat.sv_slice_vel] = vel_world_to_body_eul(Xee[cls.sv_slice_vel], Xee[cls.sv_slice_eul])
        Xbq[SixDOFBodyQuat.sv_slice_quat] = p3_alg.quat_of_euler(Xee[cls.sv_slice_eul])
        Xbq[SixDOFBodyQuat.sv_slice_rvel] = Xee[cls.sv_slice_rvel]
        return Xbq

    @classmethod
    def to_six_dof_aero_euler(cls, Xee, atm=None, t=0):
        Xae = np.zeros(SixDOFAeroEuler.sv_size)
        Xae[SixDOFAeroEuler.sv_slice_pos]   = Xee[cls.sv_slice_pos]
        wind_ned = (atm.get_wind(Xee[cls.sv_slice_pos], t) if atm is not None else [0, 0, 0])
        Xae[SixDOFAeroEuler.sv_slice_vaero] = vel_world_to_aero_eul(Xee[cls.sv_slice_vel], Xee[cls.sv_slice_eul], wind_ned)
        Xae[SixDOFAeroEuler.sv_slice_eul]   = Xee[cls.sv_slice_eul]
        Xae[SixDOFAeroEuler.sv_slice_rvel]  = Xee[cls.sv_slice_rvel]
        return Xae

    @classmethod
    def from_six_dof_aero_euler(cls, Xae, atm=None, t=0):
        return SixDOFAeroEuler.to_six_dof_euclidian_euler(Xae)

    @classmethod
    def state_from_ee(cls, Xee, atm=None, t=0.):
        return Xee
    
    @classmethod
    def cont_dyn(cls, X, t, U, P, add_weight=False):
        Fb, Mb = U[:3], U[3:]
        Xd = np.zeros(cls.sv_size)
        p_w, v_w, e_w2b, om_b = X[cls.sv_slice_pos], X[cls.sv_slice_vel], X[cls.sv_slice_eul], X[cls.sv_slice_rvel]
        # Translational kinematics
        Xd[cls.sv_slice_pos] = v_w
        # Newton for forces
        R_b2w = p3_alg.rmat_of_euler(e_w2b).T
        Xd[cls.sv_slice_vel] = 1./P.m*(np.dot(R_b2w, Fb))
        # do we add weight or not???
        if add_weight: Xd[cls.sv_slice_vel] += 1./P.m*(np.array([0, 0, P.m*P.g]))
        # Rotational kinematics
        Xd[cls.sv_slice_eul] = p3_alg.euler_derivatives(e_w2b, om_b)
        # Newton for moments
        Xd[cls.sv_slice_rvel] = np.dot(P.invI, Mb - np.cross(om_b, np.dot(P.I, om_b)))
        return Xd


    @classmethod
    def plot_trajectory(cls, time, X, U,
                        figure=None, window_title="Trajectory",
                        legend=None, label='', filename=None, atm=None):
        margins=(0.04, 0.05, 0.98, 0.96, 0.20, 0.34)
        figure = p3_pu.prepare_fig(figure, window_title, figsize=(20.48, 10.24), margins=margins)
        plots = [("x", "m", X[:,cls.sv_x]), ("y", "m", X[:,cls.sv_y]), ("z", "m", X[:,cls.sv_z]),
                 ("$\dot{x}$",  "m/s", X[:, cls.sv_xd]),
                 ("$\dot{y}$",  "m/s", X[:, cls.sv_yd]),
                 ("$\dot{z}$",  "m/s", X[:, cls.sv_zd]),
                 ("$\phi$",     "deg",   np.rad2deg(X[:,cls.sv_phi])),
                 ("$\\theta$",  "deg",   np.rad2deg(X[:,cls.sv_theta])),
                 ("$\\psi$",    "deg",   np.rad2deg(X[:,cls.sv_psi])),
                 ("$p$",        "deg/s", np.rad2deg(X[:,cls.sv_p])),
                 ("$q$",        "deg/s", np.rad2deg(X[:,cls.sv_q])),
                 ("$r$",        "deg/s", np.rad2deg(X[:,cls.sv_r]))]
        nrow = 5 if U is not None else 4
        for i, (title, ylab, data) in enumerate(plots):
            ax = plt.subplot(nrow, 3, i+1)
            plt.plot(time, data, label=label)
            p3_pu.decorate(ax, title=title, ylab=ylab, legend=legend)
        for i, min_yspan in enumerate([.5, .5, .5,  .1, .1, 1.,  1., 1., 1.,  1., 1., 1.]):
            p3_pu.ensure_yspan(plt.subplot(nrow,3,i+1), min_yspan)
        return figure

    @classmethod
    def plot_trajectory_as_ee(cls, time, Xee, U,
                              figure=None, window_title="Trajectory",
                              legend=None, label='', filename=None, atm=None):
        return cls.plot_trajectory(time, Xee, U, figure, window_title, legend, label, filename, atm)
    @classmethod
    def plot_trajectory_as_be(cls, time, Xee, U,
                              figure=None, window_title="Trajectory",
                              legend=None, label='', filename=None, atm=None):
        Xbe = np.array([cls.to_six_dof_body_euler(_X, atm, _t) for _X, _t in zip(Xee, time)])
        return SixDOFBodyEuler.plot_trajectory(time, Xbe, U, figure, window_title, legend, label, filename, atm)
    
#
#
#   AERO/Euler
#
#
class SixDOFAeroEuler():
    sv_x, sv_y, sv_z, sv_va, sv_alpha, sv_beta, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)
    sv_slice_pos   = slice(sv_x,   sv_z+1)
    sv_slice_vaero = slice(sv_va,  sv_beta+1)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)
    sv_slice_rvel  = slice(sv_p,   sv_r+1)

    @classmethod
    def get_name(cls): return 'aero_euler'

    @classmethod
    def get_default_state(cls): return np.array([0.,0.,0.,  10.,0.,0.,  0.,0.,0.,  0.,0.,0.])

    @classmethod
    def to_six_dof_euclidian_euler(cls, X, atm=None, t=0):
        Xee = np.zeros(SixDOFEuclidianEuler.sv_size)
        Xee[SixDOFEuclidianEuler.sv_slice_pos] = X[cls.sv_slice_pos]
        Xee[SixDOFEuclidianEuler.sv_slice_eul] = X[cls.sv_slice_eul]
        Xee[SixDOFEuclidianEuler.sv_slice_rvel] = X[cls.sv_slice_rvel]
        wind_ned = (atm.get_wind(X[cls.sv_slice_pos], t) if atm is not None else [0, 0, 0])
        Xee[SixDOFEuclidianEuler.sv_slice_vel] = vel_aero_to_world_euler(X[cls.sv_slice_vaero], X[cls.sv_slice_eul], wind_ned)
        return Xee

    @classmethod
    def to_six_dof_euclidian_quat(cls, Xae, atm=None, t=0):
        Xeq = np.zeros(SixDOFEuclidianQuat.sv_size)
        Xeq[SixDOFEuclidianQuat.sv_slice_pos] = Xae[cls.sv_slice_pos]
        wind_ned = (atm.get_wind(Xae[cls.sv_slice_pos], t) if atm is not None else [0, 0, 0])
        Xeq[SixDOFEuclidianEuler.sv_slice_vel] = vel_aero_to_world_euler(Xae[cls.sv_slice_vaero], Xae[cls.sv_slice_eul], wind_ned)
        Xeq[SixDOFEuclidianQuat.sv_slice_quat] = p3_alg.quat_of_euler(Xae[cls.sv_slice_eul])
        Xeq[SixDOFEuclidianEuler.sv_slice_rvel] = Xae[cls.sv_slice_rvel]
        return Xeq

    @classmethod
    def from_six_dof_aero_euler(cls, X): return X
    @classmethod
    def state_from_ee(cls, Xee, atm=None, t=0.):
        return SixDOFEuclidianEuler.to_six_dof_aero_euler(Xee, atm, t)

    @classmethod
    def state_str(cls, X):
        return """pos: {:-.2f}, {:-.2f}, {:-.2f} m
        vel: {:-.2f} m/s, alpha {:-.2f}, beta {:-.2f} deg
        att:    {:-.2f}, {:-.2f}, {:-.2f} deg
        """.format(X[cls.sv_x], X[cls.sv_y], X[cls.sv_z],
                   X[cls.sv_va], np.rad2deg(X[cls.sv_alpha]), np.rad2deg(X[cls.sv_beta]),
                   np.rad2deg(X[cls.sv_phi]), np.rad2deg(X[cls.sv_theta]), np.rad2deg(X[cls.sv_psi]))

    @classmethod
    def cont_dyn(cls, X, t, U, P, atm=None, add_weight=False):
        Fb, Mb = U[:3], U[3:]
        X_pos = X[cls.sv_slice_pos]                      # ned pos
        X_avel = va, alpha, beta = X[cls.sv_slice_vaero] # airvel, alpha, beta
        X_euler = X[cls.sv_slice_eul]                    # euler angles
        X_rvel_body = X[cls.sv_slice_rvel]               # body rotational velocities

        earth_to_body_R = p3_alg.rmat_of_euler(X_euler)
        p_w, v_aero, e_w2b, om_b = X[cls.sv_slice_pos], X[cls.sv_slice_vaero], X[cls.sv_slice_eul], X[cls.sv_slice_rvel]
        wind_ned = atm.get_wind(X_pos, t) if atm is not None else [0, 0, 0]
        waccel_body = [0, 0, 0]  # np.dot(earth_to_body, waccel_earth)
        ivel_world = vel_aero_to_world_euler(X_avel, X_euler, wind_ned)
        ivel_body = np.dot(earth_to_body_R, ivel_world)
        avel_aero = [X[cls.sv_va], 0., 0.]
        _R_aero_to_body = R_aero_to_body(alpha, beta)
        avel_body = np.dot(_R_aero_to_body, avel_aero)
        
        Xd = np.zeros(cls.sv_size)
        # Translational kinematics
        Xd[cls.sv_slice_pos] = ivel_world
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


    @classmethod
    def plot_trajectory(cls, time, X, U,
                        figure=None, window_title="Trajectory",
                        legend=None, label='', filename=None, atm=None):
        margins=(0.04, 0.05, 0.98, 0.96, 0.20, 0.34)
        figure = p3_pu.prepare_fig(figure, window_title, figsize=(20.48, 10.24), margins=margins)

        plots = [("x", "m", X[:,cls.sv_x]), ("y", "m", X[:,cls.sv_y]), ("z", "m", X[:,cls.sv_z]),
                 ("v",         "m/s", X[:,cls.sv_va]),
                 ("$\\alpha$", "deg", np.rad2deg(X[:,cls.sv_alpha])),
                 ("$\\beta$",  "deg", np.rad2deg(X[:,cls.sv_beta])),
                 ("$\phi$",    "deg", np.rad2deg(X[:,cls.sv_phi])),
                 ("$\\theta$", "deg", np.rad2deg(X[:,cls.sv_theta])),
                 ("$\\psi$",   "deg", np.rad2deg(X[:,cls.sv_psi])),
                 ("$p$",     "deg/s", np.rad2deg(X[:,cls.sv_p])),
                 ("$q$",     "deg/s", np.rad2deg(X[:,cls.sv_q])),
                 ("$r$",     "deg/s", np.rad2deg(X[:,cls.sv_r]))]
        
        nrow = 5 if U is not None else 4
        for i, (title, ylab, data) in enumerate(plots):
            ax = plt.subplot(nrow, 3, i+1)
            plt.plot(time, data, label=label)
            p3_pu.decorate(ax, title=title, ylab=ylab)
        return figure

    @classmethod
    def plot_trajectory_as_ee(cls, time, Xeq, U,
                              figure=None, window_title="Trajectory",
                              legend=None, label='', filename=None, atm=None):
        Xee = np.array([cls.to_six_dof_euclidian_euler(_X, atm, _t) for _X, _t in zip(Xeq, time)])
        return SixDOFEuclidianEuler.plot_trajectory(time, Xee, U, figure, window_title, legend, label, filename, atm)

#
#
#   Body/Euler
#
#
class SixDOFBodyEuler():
    
    sv_x, sv_y, sv_z, sv_u, sv_v, sv_w, sv_phi, sv_theta, sv_psi, sv_p, sv_q, sv_r, sv_size = range(0,13)

    sv_slice_pos   = slice(sv_x,   sv_z+1)    # position ned (north east down)
    sv_slice_vel   = slice(sv_u,   sv_w+1)    # inertial velocity in body frame (front right down)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)  # euler angles
    sv_slice_rvel  = slice(sv_p,   sv_r+1)    # rotational velocity in body frame
    @classmethod
    def get_name(cls): return 'body_euler'
    @classmethod
    def get_default_state(cls): return np.array([0.,0.,0.,  0.,0.,0.,  0.,0.,0.,  0.,0.,0.])
    @classmethod
    def to_six_dof_euclidian_euler(cls, Xbe, atm=None, t=0):
        Xee = np.zeros(SixDOFEuclidianEuler.sv_size)
        Xee[SixDOFEuclidianEuler.sv_slice_pos]  = Xbe[cls.sv_slice_pos]
        Xee[SixDOFEuclidianEuler.sv_slice_vel]  = vel_body_to_world_euler(Xbe[cls.sv_slice_vel], Xbe[cls.sv_slice_eul])
        Xee[SixDOFEuclidianEuler.sv_slice_eul]  = Xbe[cls.sv_slice_eul]
        Xee[SixDOFEuclidianEuler.sv_slice_rvel] = Xbe[cls.sv_slice_rvel]
        return Xee
    @classmethod
    def state_from_ee(cls, Xee, atm=None, t=0.):
        return SixDOFEuclidianEuler.to_six_dof_body_euler(Xee, atm, t)

    
    @classmethod
    def cont_dyn(cls, X, t, U, P, add_weight=False):
        Fb, Mb = U[:3], U[3:]
        Xd = np.zeros(cls.sv_size)
        p_w, ivel_b, e_w2b, rvel_b = X[cls.sv_slice_pos], X[cls.sv_slice_vel], X[cls.sv_slice_eul], X[cls.sv_slice_rvel]
        R_b2w = p3_alg.rmat_of_euler(e_w2b).T
        # Translational kinematics
        Xd[cls.sv_slice_pos] = vel_body_to_world_R(ivel_b, R_b2w)
        # Newton for forces  # FIXME  why 2*np.cross?
        Xd[cls.sv_slice_vel] = 1/P.m*Fb-2.*np.cross(rvel_b, ivel_b) ### CHECK with silvia
        # do we add weight or not???
        #if add_weight: Xd[cls.sv_slice_vel] += 1./P.m*(np.array([0, 0, P.m*P.g]))
        # Rotational kinematics
        Xd[cls.sv_slice_eul] = p3_alg.euler_derivatives(e_w2b, rvel_b)
        # Newton for moments
        Xd[cls.sv_slice_rvel] = np.dot(P.invI, Mb - np.cross(rvel_b, np.dot(P.I, rvel_b)))
        return Xd

    @classmethod
    def plot_trajectory(cls, time, X, U,
                        figure=None, window_title="Trajectory be",
                        legend=None, label='', filename=None, atm=None):
        margins=(0.04, 0.05, 0.98, 0.96, 0.20, 0.34)
        figure = p3_pu.prepare_fig(figure, window_title, figsize=(20.48, 10.24), margins=margins)
        plots = [("x", "m", X[:,cls.sv_x]), ("y", "m", X[:,cls.sv_y]), ("z", "m", X[:,cls.sv_z]),
                 ("$u$",        "m/s", X[:, cls.sv_u]),
                 ("$v$",        "m/s", X[:, cls.sv_v]),
                 ("$w$",        "m/s", X[:, cls.sv_w]),
                 ("$\phi$",     "deg",   np.rad2deg(X[:,cls.sv_phi])),
                 ("$\\theta$",  "deg",   np.rad2deg(X[:,cls.sv_theta])),
                 ("$\\psi$",    "deg",   np.rad2deg(X[:,cls.sv_psi])),
                 ("$p$",        "deg/s", np.rad2deg(X[:,cls.sv_p])),
                 ("$q$",        "deg/s", np.rad2deg(X[:,cls.sv_q])),
                 ("$r$",        "deg/s", np.rad2deg(X[:,cls.sv_r]))]
        nrow = 5 if U is not None else 4
        for i, (title, ylab, data) in enumerate(plots):
            ax = plt.subplot(nrow, 3, i+1)
            plt.plot(time, data, label=label)
            p3_pu.decorate(ax, title=title, ylab=ylab, legend=legend)
        for i, min_yspan in enumerate([.5, .5, .5,  .1, .1, 1.,  1., 1., 1.,  1., 1., 1.]):
            p3_pu.ensure_yspan(plt.subplot(nrow,3,i+1), min_yspan)
        return figure
    @classmethod
    def plot_trajectory_as_be(cls, time, Xbe, U,
                              figure=None, window_title="Trajectory",
                              legend=None, label='', filename=None, atm=None):
        return cls.plot_trajectory(time, Xbe, U, figure, window_title, legend, label, filename, atm)
    @classmethod
    def plot_trajectory_as_ee(cls, time, Xbe, U,
                              figure=None, window_title="Trajectory",
                              legend=None, label='', filename=None, atm=None):
        Xee = np.array([cls.to_six_dof_euclidian_euler(_X, atm, _t) for _X, _t in zip(Xbe, time)])
        return SixDOFEuclidianEuler.plot_trajectory(time, Xee, U, figure, window_title, legend, label, filename, atm)

    

#
#
#   Euclidian/Quat
#
#
class SixDOFEuclidianQuat():
    
    sv_x, sv_y, sv_z, sv_xd, sv_yd, sv_zd, sv_qi, sv_qx, sv_qy, sv_qz, sv_p, sv_q, sv_r, sv_size = range(0,14)
    
    sv_slice_pos   = slice(sv_x,  sv_z+1)   # position ned (north east down)
    sv_slice_vel   = slice(sv_xd, sv_zd+1)  # inertial velocity in world frame (ned)
    sv_slice_quat  = slice(sv_qi, sv_qz+1)  # quaternion from world to body ( FIXME check direction )
    sv_slice_rvel  = slice(sv_p,  sv_r+1)   # rotational velocity in body frame

    @classmethod
    def get_name(cls): return 'euclidian_quat'
    @classmethod
    def get_default_state(cls): return np.array([0.,0.,0.,  0.,0.,0.,  1.,0.,0.,0.,  0.,0.,0.])
    
    @classmethod
    def from_six_dof_aero_euler(cls, Xae, atm=None, t=0):
        return SixDOFAeroEuler.to_six_dof_euclidian_quat(Xae, atm, t)
    @classmethod
    def state_from_ee(cls, Xee, atm=None, t=0.):
        return SixDOFEuclidianEuler.to_six_dof_euclidian_quat(Xee, atm, t)

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

    @classmethod
    def plot_trajectory(cls, time, Xeq, U,
                        figure=None, window_title="Trajectory",
                        legend=None, label='', filename=None, atm=None):
        return cls.plot_trajectory_as_ee(time, Xeq, U, figure, window_title, legend, label, filename, atm)

    @classmethod
    def plot_trajectory_as_ee(cls, time, Xeq, U,
                              figure=None, window_title="Trajectory",
                              legend=None, label='', filename=None, atm=None):
        Xee = np.array([cls.to_six_dof_euclidian_euler(_X, atm, _t) for _X, _t in zip(Xeq, time)])
        return SixDOFEuclidianEuler.plot_trajectory(time, Xee, U, figure, window_title, legend, label, filename, atm)
#
#
#   Body/Quat
#
#
class SixDOFBodyQuat():
    
    sv_x, sv_y, sv_z, sv_u, sv_v, sv_w, sv_qi, sv_qx, sv_qy, sv_qz, sv_p, sv_q, sv_r, sv_size = range(0,14)
    
    sv_slice_pos   = slice(sv_x,  sv_z+1)   # position ned (north east down)
    sv_slice_vel   = slice(sv_u, sv_w+1)    # inertial velocity in body frame (front right down)
    sv_slice_quat  = slice(sv_qi, sv_qz+1)  # quaternion from world to body ( FIXME check direction )
    sv_slice_rvel  = slice(sv_p,  sv_r+1)   # rotational velocity in body frame

    @classmethod
    def get_name(cls): return 'body_quat'
    @classmethod
    def get_default_state(cls): return np.array([0.,0.,0.,  0.,0.,0.,  1.,0.,0.,0.,  0.,0.,0.])
    
    
    @classmethod
    def to_six_dof_euclidian_euler(cls, Xbq, atm=None, t=0):
        Xee = np.zeros(SixDOFEuclidianEuler.sv_size)
        Xee[SixDOFEuclidianEuler.sv_slice_pos]  = Xbq[cls.sv_slice_pos]
        Xee[SixDOFEuclidianEuler.sv_slice_vel]  = vel_body_to_world_quat(Xbq[cls.sv_slice_vel], Xbq[cls.sv_slice_quat])
        Xee[SixDOFEuclidianEuler.sv_slice_eul]  = p3_alg.euler_of_quat(Xbq[cls.sv_slice_quat])
        Xee[SixDOFEuclidianEuler.sv_slice_rvel] = Xbq[cls.sv_slice_rvel]
        return Xee
    @classmethod
    def to_six_dof_body_euler(cls, Xbq, atm=None, t=0):
        Xbe = np.zeros(SixDOFBodyEuler.sv_size)
        Xbe[SixDOFEuclidianEuler.sv_slice_pos]  = Xbq[cls.sv_slice_pos]
        Xbe[SixDOFEuclidianEuler.sv_slice_vel]  = Xbq[cls.sv_slice_vel]
        Xbe[SixDOFEuclidianEuler.sv_slice_eul]  = p3_alg.euler_of_quat(Xbq[cls.sv_slice_quat])
        Xbe[SixDOFEuclidianEuler.sv_slice_rvel] = Xbq[cls.sv_slice_rvel]
        return Xbe
    @classmethod
    def state_from_ee(cls, Xee, atm=None, t=0.):
        return SixDOFEuclidianEuler.to_six_dof_body_quat(Xee, atm, t)

    @classmethod
    def cont_dyn(cls, X, t, U, P, add_weight=False):
        Fb, Mb = U[:3], U[3:]
        Xd = np.zeros(cls.sv_size)
        p_w, ivel_b, q_w2b, rvel_b = X[cls.sv_slice_pos], X[cls.sv_slice_vel], X[cls.sv_slice_quat], X[cls.sv_slice_rvel]
        R_b2w = p3_alg.rmat_of_quat(q_w2b).T
        # Translational kinematics
        Xd[cls.sv_slice_pos] = vel_body_to_world_R(ivel_b, R_b2w)
        # Newton for forces
        Xd[cls.sv_slice_vel] = 1/P.m*Fb-2.*np.cross(rvel_b, ivel_b)
        # do we add weight or not???
        #if add_weight: Xd[cls.sv_slice_vel] += 1./P.m*(np.array([0, 0, P.m*P.g]))
        # Rotational kinematics
        Xd[cls.sv_slice_quat] = p3_alg.quat_derivative(q_w2b, rvel_b)
        # Newton for moments
        Xd[cls.sv_slice_rvel] = np.dot(P.invI, Mb - np.cross(rvel_b, np.dot(P.I, rvel_b)))
        return Xd

    @classmethod
    def plot_trajectory(cls, time, Xbq, U,
                        figure=None, window_title="Trajectory",
                        legend=None, label='', filename=None, atm=None):
        # FIXME, plot uvw
        return cls.plot_trajectory_as_ee(time, Xbq, U, figure, window_title, legend, label, filename, atm)
    @classmethod
    def plot_trajectory_as_be(cls, time, Xbq, U,
                              figure=None, window_title="Trajectory",
                              legend=None, label='', filename=None, atm=None):
        Xbe = np.array([cls.to_six_dof_body_euler(_X, atm, _t) for _X, _t in zip(Xbq, time)])
        return SixDOFBodyEuler.plot_trajectory(time, Xbe, U, figure, window_title, legend, label, filename, atm)
    @classmethod
    def plot_trajectory_as_ee(cls, time, Xbq, U,
                              figure=None, window_title="Trajectory",
                              legend=None, label='', filename=None, atm=None):
        Xee = np.array([cls.to_six_dof_euclidian_euler(_X, atm, _t) for _X, _t in zip(Xbq, time)])
        return SixDOFEuclidianEuler.plot_trajectory(time, Xee, U, figure, window_title, legend, label, filename, atm)
