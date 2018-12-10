import numpy as np


class Sim:

    def __init__(self, fdm, ctl):
        self.dt_internal = 0.01
        self.fdm, self.ctl = fdm, ctl
        self.Yc = np.zeros(ctl.spv_size)
        
    def reset(self, t):
        self.t = t
        self.Xe, self.Ue = self.fdm.trim()
        self.X0 = np.array(self.Xe)
        #self.X0[fdm.sv_y] += 0.1
        #self.X0[fdm.sv_zd] += 0.05
        #self.X0[fdm.sv_r] += 0.05
        #phi0, theta0, psi0 = np.deg2rad(10.), 0, 0
        #self.X0[fdm.sv_slice_quat] = pal.quat_of_euler([phi0, theta0, psi0])
        return self.fdm.reset(self.X0, self.t)
        
    def run(self, t):
        U = self.ctl.get(t, self.fdm.X, self.Yc)
        while t - self.t > 1e-6:
            U = self.ctl.get(t, self.fdm.X, self.Yc)
            self.t += self.dt_internal
            self.fdm.run(self.t, U)
        return U, self.fdm.X
        
