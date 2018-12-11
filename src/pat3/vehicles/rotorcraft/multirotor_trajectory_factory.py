import numpy as np

import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt



def traj1(): return pmt.Line([0, 0, 0], [1, 1, -0.5], v=0.1)

def traj2(): return pmt.Circle([0, 0, -0.5], r=1., v=1.)

def traj3(): return pmt.Circle([0, 0, -0.5], r=1., v=4.)

def traj4(): return pmt.Circle([0, 0, -0.5], r=2., v=1.)

def traj5():
    r, v = 1., 4.; om = v/r
    return pmt.Circle(c=[0, 0, -0.5], r=1., v=4., alpha0=0., dalpha=2*np.pi, psit=pmt.SinOne(om=om))

def traj6(): return pmt.Oval(l=2, r=1, v=4)

def traj7(): return pmt.DoubleOval(l=2, r=1, v=4)

def traj8(): return pmt.FigureOfEight()

trajectories = {'line1':   ('line v=0.1',      traj1),
                'circle1': ('circle r=1 v=1',  traj2),
                'circle2': ('circle r=1 v=4',  traj3),
                'circle3': ('circle r=2 v=4',  traj4),
                'circle4': ('circle r=1 v=4 looking at center',  traj5),
                'oval1':   ('oval',            traj6),
                'doval1':  ('double oval',     traj7),
                'foe1':    ('figure of eight', traj8)
}


def get(traj_name):
    return trajectories[traj_name][1](), trajectories[traj_name][0]

def list():
    names = ['{}: {}'.format(k,v[0]) for k,v in sorted(trajectories.iteritems())]
    return names
