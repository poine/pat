import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

import pat3.algebra    as pal
import pat3.plot_utils as ppu
import pat3.flat_world as pfw


class DynamicModel:

    iv_size = 1
    sv_size = 1
    
    def __init__(self, sv_size, P=None, dt=0.005, solve_ode_first_order=False):
        self.sv_size = sv_size
        self.dt, self.P = dt, P
        self.solve_ode_first_order = solve_ode_first_order
        self.t, self.X = 0., np.zeros((self.sv_size))
        # byproducts
        self.T_w2b = np.eye(4)  # world (ned) to body (frd) homogeneous transform

    def disc_dyn(self, Xk, tk, Uk, dt):
        '''
        Discrete-time State Space Representation: Xk+1 = f_param(Xk, Uk)
        '''
        if self.solve_ode_first_order:
            Xkp1 = Xk + self.cont_dyn(Xk, tk, Uk)*dt # first order
        else:
            _unused, Xkp1 = scipy.integrate.odeint(self.cont_dyn, Xk, [tk, tk+dt], args=(Uk, ))
            
        return Xkp1

    # def reset(self, X0, t0):
    #     self.X, self.t = X0, t0
    #     self.update_byproducts()
    #     return self.X
    
    # def run(self, tf, U):
    #     remaining_to_tf = tf - self.t
    #     while remaining_to_tf > 0:
    #         dt = min(self.dt, remaining_to_tf)
    #         #print(' integrate fdm from {} to {}'.format(self.t, self.t+dt))
    #         self.X = self.disc_dyn(self.X, self.t, U, dt)
    #         self.t += dt
    #         remaining_to_tf = tf - self.t
    #     self.update_byproducts()
    #     return self.X
    
    # def update_byproducts(self):
    #     self.T_w2b[:3,:3] = pal.rmat_of_quat(self.X[sv_slice_quat]).T # that is freaking weird....
    #     self.T_w2b[:3,3] = self.X[sv_slice_pos]


        
class SolidParam:
    def __init__(self):
        self.m = 1.      # mass
        self.g = 9.81    # gravity FIXME, should not be constant
        self.J = np.diag([0.01, 0.01, 0.05]) # inertia
        # precompute
        self.invJ = np.linalg.inv(self.J)  # inertia invert
        
class SolidFDM(DynamicModel):
    ''' 
    An object encapsulating the dynamic model of a solid
    '''

    '''State vector components'''
    sv_x    = 0  # position in euclidian world
    sv_y    = 1  #
    sv_z    = 2  # velocity in euclidian world
    sv_xd   = 3  #
    sv_yd   = 4  #
    sv_zd   = 5  #
    sv_qi   = 6  # world to body rotation quaternion
    sv_qx   = 7  #
    sv_qy   = 8  #
    sv_qz   = 9  #
    sv_p   = 10  # rotational velocities in body frame
    sv_q   = 11  # ( front, right, down )
    sv_r   = 12  #
    sv_size = 13

    sv_slice_pos  = slice(sv_x, sv_z+1)
    sv_slice_vel  = slice(sv_xd, sv_zd+1)
    sv_slice_quat = slice(sv_qi, sv_qz+1)
    sv_slice_rvel = slice(sv_p,sv_r+1)

    iv_fxb   = 0  # forces in body frame
    iv_fyb   = 1
    iv_fzb   = 2
    iv_mxb   = 3  # moments in body frame
    iv_myb   = 4
    iv_mzb   = 5
    iv_size =  6  # size of input vector
    
    def __init__(self, P=None, dt=0.005):
        P = P if P is not None else SolidParam()
        DynamicModel.__init__(self, self.sv_size, P, dt, solve_ode_first_order=False)
        
    def cont_dyn(self, X, t, U):
        Fb, Mb = U[:3], U[3:]
        Xd = np.zeros(self.sv_size)
        p_w, v_w, q_w2b, om_b = X[self.sv_slice_pos], X[self.sv_slice_vel], X[self.sv_slice_quat], X[self.sv_slice_rvel]
        # Translational kinematics
        Xd[self.sv_slice_pos] = v_w
        # Newton for forces
        R_w2b =  pal.rmat_of_quat(q_w2b)
        Xd[self.sv_slice_vel] = 1./self.P.m*(np.dot(R_w2b.T, Fb) + [0, 0, self.P.m*self.P.g])
        #Xd[self.sv_slice_vel] = 1./self.P.m*(np.dot(R_w2b.T, Fb))
        # Rotational kinematics
        Xd[self.sv_slice_quat] = pal.quat_derivative(q_w2b, om_b)
        # Newton for moments
        Xd[self.sv_slice_rvel] = np.dot(self.P.invJ, Mb - np.cross(om_b, np.dot(self.P.J, om_b)))
        return Xd

    def trim(self):
        Xe = np.zeros(self.sv_size); Xe[self.sv_qi] = 1.
        Fb, Mb = [0, 0, -self.P.m*self.P.g], [0, 0, 0]
        Ue = np.concatenate((Fb, Mb))
        return Xe, Ue

    def plot(self, time, X, U=None, figure=None, window_title="Solid Trajectory"):
        figure = ppu.prepare_fig(figure, window_title, (20.48, 10.24))
 
        eulers = np.array([pal.euler_of_quat(_q) for _q in X[:,self.sv_slice_quat]])
        phi, theta, psi = eulers[:,0], eulers[:,1], eulers[:,2]
        plots = [("$x$",       "m",     0.5, X[:,self.sv_x]),
                 ("$y$",       "m",     0.5, X[:,self.sv_y]),
                 ("$z$",       "m",     0.5, X[:,self.sv_z]),
                 ("$\dot{x}$", "m/s",   0.5, X[:,self.sv_xd]),
                 ("$\dot{y}$", "m/s",   0.5, X[:,self.sv_yd]),
                 ("$\dot{z}$", "m/s",   0.5, X[:,self.sv_zd]),
                 ("$\phi$",    "deg",   0.5, np.rad2deg(phi)),
                 ("$\\theta$", "deg",   0.5, np.rad2deg(theta)),
                 ("$\psi$",    "deg",   0.5, np.rad2deg(psi)),
                 ("$p$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_p])),
                 ("$q$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_q])),
                 ("$r$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_r])),
        ]
        if U is not None: # force an extra row of plots
            foo = np.empty((len(time))); foo.fill(np.nan)
            plots += [("$Fb$", "N", 0.1, foo), ("$Mb$", "Nm", 0.1, foo)]
        figure = ppu.plot_in_grid(time, plots, 3, figure, window_title)
        if U is not None:
            ax = plt.subplot(5, 3, 13)
            for i,txt in enumerate(('fx', 'fy', 'fz')):
                plt.plot(time, U[:,i], label=txt)
                plt.legend()
            ax = plt.subplot(5, 3, 14)
            for i,txt in enumerate(('mx', 'my', 'mz')):
                plt.plot(time, U[:,i+3], label=txt)
                plt.legend()
        return figure

#
#
# Using euler angle and body velocity: unfinished/untested
#
#
class SolidDM2(DynamicModel):
    '''
    State vector components
    '''
    sv_x    = 0  # position in euclidian world (north east down)
    sv_y    = 1  #
    sv_z    = 2  # velocity in body frame
    sv_u    = 3  # ( front, right, down )
    sv_v    = 4  #
    sv_w    = 5  #
    sv_phi  = 6  # euler angles
    sv_the  = 7  #
    sv_psi  = 8  # 
    sv_p    = 9  # rotational velocities in body frame
    sv_q    = 10 # ( front, right, down )
    sv_r    = 11 #
    sv_size = 12
    
    sv_slice_pos  = slice(sv_x, sv_z+1)
    sv_slice_uvw  = slice(sv_u, sv_w+1)
    sv_slice_eul  = slice(sv_phi, sv_psi+1)
    sv_slice_rvel = slice(sv_p,sv_r+1)

    iv_fxb   = 0  # forces in body frame
    iv_fyb   = 1
    iv_fzb   = 2
    iv_mxb   = 3  # moments in body frame
    iv_myb   = 4
    iv_mzb   = 5
    iv_size =  6  # size of input vector

    def __init__(self, P=None, dt=0.005):
        P = P if P is not None else SolidParam()
        DynamicModel.__init__(self, self.sv_size, P, dt, solve_ode_first_order=False)
        # byproducts
        self.T_w2b = np.eye(4)  # world (ned) to body (frd) homogeneous transform

    def trim(self):
        Xe = np.zeros(self.sv_size)
        Fb, Mb = [0, 0, -self.P.m*self.P.g], [0, 0, 0]
        Ue = np.concatenate((Fb, Mb))
        return Xe, Ue

    def cont_dyn(self, X, t, U):
        Fb, Mb = U[:3], U[3:]
        vel_body, euler, rvel_body = X[self.sv_slice_uvw], X[self.sv_slice_eul], X[self.sv_slice_rvel]
        body_to_earth = pfw.body_to_earth(euler)
        Xdot = np.zeros(self.sv_size)
        Xdot[self.sv_slice_pos] = np.dot(body_to_earth, vel_body)
        inv_mamat = None
        if inv_mamat != None:
            Xdot[s1_u:s1_w+1] = np.dot(inv_mamat, Fb -2.*np.cross(rvel_body,  vel_body))
        else:
            Xdot[self.sv_slice_uvw] = 1/self.P.m*Fb-2.*np.cross(rvel_body,  vel_body)
        Xdot[self.sv_slice_eul] = pal.euler_derivatives(euler, rvel_body)
        raccel_body = np.dot(self.P.invJ, Mb - np.cross(rvel_body, np.dot(self.P.J, rvel_body)))
        Xdot[self.sv_slice_rvel] = raccel_body
        return Xdot
        # where does that come from???
        # Fb, Mb = U[:3], U[3:]
        # Xd = np.zeros(self.sv_size)
        # p_w, v_w, q_w2b, om_b = X[self.sv_slice_pos], X[self.sv_slice_vel], X[self.sv_slice_quat], X[self.sv_slice_rvel]
        # # Translational kinematics
        # Xd[self.sv_slice_pos] = v_w
        # # Newton for forces
        # R_w2b =  pal.rmat_of_quat(q_w2b)
        # Xd[self.sv_slice_vel] = 1./self.P.m*(np.dot(R_w2b.T, Fb) + [0, 0, self.P.m*self.P.g])
        # # Rotational kinematics
        # Xd[self.sv_slice_quat] = pal.quat_derivative(q_w2b, om_b)
        # # Newton for moments
        # Xd[self.sv_slice_rvel] = np.dot(self.P.invJ, Mb - np.cross(om_b, np.dot(self.P.J, om_b)))
        # return Xd

    def plot(self, time, X, U=None, figure=None, window_title="Trajectory"):
        figure = ppu.prepare_fig(figure, window_title, (20.48, 10.24))
        phi, theta, psi = X[:,self.sv_phi], X[:,self.sv_the], X[:,self.sv_psi]
        plots = [("$x$",       "m",     0.5, X[:,self.sv_x]),
                 ("$y$",       "m",     0.5, X[:,self.sv_y]),
                 ("$z$",       "m",     0.5, X[:,self.sv_z]),
                 ("$u$",       "m/s",   0.5, X[:,self.sv_u]),
                 ("$v$",       "m/s",   0.5, X[:,self.sv_v]),
                 ("$w$",       "m/s",   0.5, X[:,self.sv_w]),
                 ("$\phi$",    "deg",   0.5, np.rad2deg(phi)),
                 ("$\\theta$", "deg",   0.5, np.rad2deg(theta)),
                 ("$\psi$",    "deg",   0.5, np.rad2deg(psi)),
                 ("$p$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_p])),
                 ("$q$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_q])),
                 ("$r$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_r])),
        ]
        if U is not None: # force an extra row of plots
            foo = np.empty((len(time))); foo.fill(np.nan)
            plots += [("$Fb$", "N", 0.1, foo), ("$Mb$", "Nm", 0.1, foo)]
        figure = ppu.plot_in_grid(time, plots, 3, figure, window_title)
        if U is not None:
            ax = plt.subplot(5, 3, 13)
            for i,txt in enumerate(('fx', 'fy', 'fz')):
                plt.plot(time, U[:,i], label=txt)
                plt.legend()
            ax = plt.subplot(5, 3, 14)
            for i,txt in enumerate(('mx', 'my', 'mz')):
                plt.plot(time, U[:,i+3], label=txt)
                plt.legend()
        return figure



#
#
# Using euler angles and aero velocity
#
#
class SolidDM3(DynamicModel):
    
    """
    Components of the state vector
    """
    sv_x     = 0   # position x axis
    sv_y     = 1   # position y axis
    sv_z     = 2   # heigh above ground
    sv_v     = 3   # airspeed
    sv_alpha = 4   # alpha
    sv_beta  = 5   # beta
    sv_phi   = 6   # roll  (euler, ltp to body)
    sv_theta = 7   # pitch (euler, ltp to body)
    sv_psi   = 8   # yaw   (euler, ltp to body)
    sv_p     = 9   # rotational vel body x
    sv_q     = 10  # rotational vel body y
    sv_r     = 11  # rotational vel body z
    sv_size  = 12
    
    sv_slice_pos   = slice(sv_x, sv_z+1)
    sv_slice_vaero = slice(sv_v, sv_beta+1)
    sv_slice_eul   = slice(sv_phi, sv_psi+1)
    sv_slice_rvel  = slice(sv_p,sv_r+1)

    iv_fxb   = 0  # forces in body frame
    iv_fyb   = 1
    iv_fzb   = 2
    iv_mxb   = 3  # moments in body frame
    iv_myb   = 4
    iv_mzb   = 5
    iv_size =  6  # size of input vector

    def __init__(self, P=None, dt=0.005):
        P = P if P is not None else SolidParam()
        DynamicModel.__init__(self, self.sv_size, P, dt, solve_ode_first_order=False)
        # byproducts
        self.T_w2b = np.eye(4)  # world (ned) to body (frd) homogeneous transform

    def trim(self):
        Xe = np.zeros(self.sv_size)
        Fb, Mb = [0, 0, -self.P.m*self.P.g], [0, 0, 0]
        Ue = np.concatenate((Fb, Mb))
        return Xe, Ue

    def cont_dyn(self, X, t, U):        
        Fb, Mb = U[:3], U[3:]

        Xdot = np.zeros(self.sv_size)
        X_euler, X_rvel = X[self.sv_slice_eul], X[self.sv_slice_rvel]
        body_to_earth_R = pfw.body_to_earth(X_euler)
        aero_to_body_R = pfw.aero_to_body(X[self.sv_alpha], X[self.sv_beta])
        vel_body = u,v,w = np.dot(aero_to_body_R, [X[self.sv_v], 0., 0.]) 
        # translational kinematics
        Xdot[self.sv_slice_pos] = np.dot(body_to_earth_R, vel_body)
        # rotational kinematics
        Xdot[self.sv_slice_eul] = pal.euler_derivatives(X_euler, X_rvel)
        # translational dynamics
        accel_body = ud, vd, wd = 1./self.P.m*Fb - np.cross(X_rvel, vel_body)
        Xdot[self.sv_v] = np.inner(vel_body, accel_body)/X[self.sv_v]
        Xdot[self.sv_alpha] = (u*wd - w*ud)/(u**2+w**2)
        Xdot[self.sv_beta] = (X[self.sv_v]*vd - v*Xdot[self.sv_v]) / X[self.sv_v] / np.sqrt(u**2+w**2)
        # rotational dynamics
        raccel_body = np.dot(self.P.invJ, Mb - np.cross(X_rvel, np.dot(self.P.J, X_rvel)))
        Xdot[self.sv_slice_rvel] = raccel_body
        return Xdot

    def plot(self, time, X, U=None, figure=None, window_title="Trajectory", extra_rows=0):
        figure = ppu.prepare_fig(figure, window_title, (20.48, 10.24))
        phi, theta, psi = X[:,self.sv_phi], X[:,self.sv_theta], X[:,self.sv_psi]
        plots = [("$x$",       "m",     0.5, X[:,self.sv_x]),
                 ("$y$",       "m",     0.5, X[:,self.sv_y]),
                 ("$z$",       "m",     0.5, X[:,self.sv_z]),
                 ("$v$",       "m/s",   0.5, X[:,self.sv_v]),
                 ("$\\alpha$",  "deg",   0.5, np.rad2deg(X[:,self.sv_alpha])),
                 ("$\\beta$",   "deg",   0.5, np.rad2deg(X[:,self.sv_beta])),
                 ("$\phi$",    "deg",   0.5, np.rad2deg(phi)),
                 ("$\\theta$", "deg",   0.5, np.rad2deg(theta)),
                 ("$\psi$",    "deg",   0.5, np.rad2deg(psi)),
                 ("$p$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_p])),
                 ("$q$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_q])),
                 ("$r$",       "deg/s", 0.5, np.rad2deg(X[:,self.sv_r])),
        ]
        # if U is not None: # force an extra row of plots
        #     foo = np.empty((len(time))); foo.fill(np.nan)
        #     plots += [("$Fb$", "N", 0.1, foo), ("$Mb$", "Nm", 0.1, foo)]
        figure = ppu.plot_in_grid(time, plots, 3, figure, window_title, extra_rows=1+extra_rows if U is not None else extra_rows)
        if U is not None:
            ax = plt.subplot(5, 3, 13)
            for i,txt in enumerate(('fx', 'fy', 'fz')):
                plt.plot(time, U[:,i], label=txt)
                plt.legend()
            ax = plt.subplot(5, 3, 14)
            for i,txt in enumerate(('mx', 'my', 'mz')):
                plt.plot(time, U[:,i+3], label=txt)
                plt.legend()
        return figure



#
# A solid that compensates gravity
#
class UfoParam(SolidParam):
    def __init__(self):
        SolidParam.__init__(self)
        self.Cd = 0.1
        
class UFOFDM(SolidFDM):

    iv_fz   = 0
    iv_mx   = 1
    iv_my   = 2
    iv_mz   = 3
    iv_size = 4
    
    def __init__(self, P=None, dt=0.005):
        P = P if P is not None else UfoParam()
        SolidFDM.__init__(self, P, dt)
    
    def cont_dyn(self, X, t, U):
        Fb, Mb = [0, 0, -U[0]], U[1:] # thrust, Moments
        Dw = -self.P.Cd*X[self.sv_slice_vel]
        R_w2b = pal.rmat_of_quat(X[self.sv_slice_quat])
        Fb += np.dot(R_w2b, Dw) # drag
        Fb -= np.dot(R_w2b, [0., 0., self.P.m*self.P.g]) # lift
        return SolidFDM.cont_dyn(self, X, t, np.concatenate((Fb, Mb)))
        
    def trim(self):
        Xe = np.zeros(self.sv_size); Xe[self.sv_qi] = 1.
        Ue = np.array([0, 0, 0, 0])
        return Xe, Ue

    def plot(self, time, X, U=None, figure=None, window_title="UFO Trajectory"):
        Usolid = np.zeros((len(time), SolidFDM.iv_size))
        Usolid[:,SolidFDM.iv_fzb] = U[:,self.iv_fz]
        Usolid[:,SolidFDM.iv_mxb] = U[:,self.iv_mx]
        Usolid[:,SolidFDM.iv_myb] = U[:,self.iv_my]
        Usolid[:,SolidFDM.iv_mzb] = U[:,self.iv_mz]
        return SolidFDM.plot(self, time, X, Usolid, figure, window_title)
