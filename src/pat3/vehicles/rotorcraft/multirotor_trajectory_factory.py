import numpy as np

import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt


class Traj1(pmt.Line):
    name, desc = 'line1', 'line v=0.1'
    def __init__(self): pmt.Line.__init__(self, [0, 0, 0], [1, 1, -0.5], v=0.1)

class Traj2(pmt.Circle):
    name, desc = 'circle1', 'circle r=1 v=1'
    def __init__(self): pmt.Circle.__init__(self, [0, 0, -0.5], r=1., v=1.)

class Traj3(pmt.Circle):
    name,desc = 'circle2', 'circle r=1 v=4'
    def __init__(self): pmt.Circle.__init__(self, [0, 0, -0.5], r=1., v=4.)

class Traj4(pmt.Circle):
    name, desc = 'circle3', 'circle r=2 v=4'
    def __init__(self): pmt.Circle.__init__(self, [0, 0, -0.5], r=2., v=1.)

class Traj5(pmt.Circle):
    name, desc = 'circle3', 'circle r=2 v=4'
    def __init__(self):
        r, v = 1., 4.; om = v/r
        pmt.Circle.__init__(self, c=[0, 0, -0.5], r=1., v=4., alpha0=0., dalpha=2*np.pi, psit=pmt.SinOne(om=om))

class Traj6(pmt.Oval):
    name, desc = 'oval1', 'oval'
    def __init__(self): pmt.Oval.__init__(self, l=2, r=1, v=4)

class Traj7(pmt.DoubleOval):
    name, desc = 'doval1', 'double oval'
    def __init__(self): return pmt.DoubleOval.__init__(self, l=2, r=1, v=4)
    
class Traj8(pmt.FigureOfEight):
    name, desc = 'foe1', 'figure of eight'

class Traj9(pmt.Foo):
    name, desc = 'foo', 'foo'

    


trajectories = {T.name: (T.desc, T) for T in [Traj1, Traj2, Traj3, Traj4, Traj5, Traj6, Traj7, Traj8, Traj9]}

def get(traj_name):
    return trajectories[traj_name][1](), trajectories[traj_name][0]

def list():
    names = ['{}: {}'.format(k,v[0]) for k,v in sorted(trajectories.iteritems())]
    return names
