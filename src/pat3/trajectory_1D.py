#
# 1D trajectories
#
import numpy as np

# (we need to move some stuff from rotorcraft trajectories in here
import pat3.vehicles.rotorcraft.multirotor_trajectory as ptraj

class AffOne:
    def __init__(self, P0, P1):
        self.P0, self.P1 = P0, P1
        self.duration = P1[0]-P0[0]
        self.c = (P1[1]-P0[1])/self.duration

    def get(self, t):
        y = self.P0[1] + (t-self.P0[0])*self.c
        return [y, self.c, 0, 0, 0]

    
class CompositeOne:
    def __init__(self, trajs):
        self.trajs = trajs
        self.durations = [trj.duration for trj in self.trajs]
        self.ends = np.cumsum(self.durations)
        self.duration = self.ends[-1]
        
    def get(self, t):
        cur_step = np.argmax(self.ends >= t)
        return self.trajs[cur_step].get(t)

class SmoothedCompositeOne(CompositeOne):
    def __init__(self, trajs, eps=0.15):
        CompositeOne.__init__(self, trajs)
        self.eps = eps
        self.corners = []
        for i in range(len(self.trajs)-1):
            Y0, Y1 = self.trajs[i].get(self.ends[i]-self.eps), self.trajs[i+1].get(self.ends[i]+self.eps)
            self.corners.append(ptraj.PolynomialOne(Y0, Y1, 2*self.eps))
        
    def get(self, t):
        cur_step = np.argmax(self.ends >= t)
        if cur_step<len(self.trajs)-1 and t>self.ends[cur_step]-self.eps:
            return self.corners[cur_step].get(t-self.ends[cur_step]+self.eps)
        elif cur_step > 0 and t<self.ends[cur_step-1]+self.eps:
            return self.corners[cur_step-1].get(t-self.ends[cur_step-1]+self.eps)
        else: return self.trajs[cur_step].get(t)

