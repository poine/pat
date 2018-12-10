import numpy as np

import pdb

class Sim:

    def __init__(self, fdm, ctl):
        self.fdm, self.ctl = fdm, ctl
        self.Yc = np.zeros(ctl.spv_size)
        
    def reset(self, t0, X0=None):
        if X0 is None:
            X0, U0 = self.fdm.trim()
        return self.fdm.reset(X0, t0)
        
    def run(self, t1):
        U = self.ctl.get(self.fdm.t, self.fdm.X, self.Yc) # in case we don't run the fdm
        #pdb.set_trace()
        #print('run sim to {}'.format(t1))
        while t1 - self.fdm.t > 0:#self.fdm.dt:
            #print(' compute control at {}'.format(self.fdm.t))
            U = self.ctl.get(self.fdm.t, self.fdm.X, self.Yc)
            self.fdm.run(self.fdm.t+self.fdm.dt, U)
        return U, self.fdm.X
        
