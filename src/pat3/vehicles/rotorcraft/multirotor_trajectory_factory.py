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
    def __init__(self): pmt.Oval.__init__(self, l=1, r=1, v=3., z=-0.5)

class Traj7(pmt.DoubleOval):
    name, desc = 'doval1', 'double oval'
    def __init__(self): pmt.DoubleOval.__init__(self, l=2, r=1, v=4)
    
class Traj8(pmt.FigureOfEight):
    name, desc = 'foe1', 'figure of eight'

class Traj9(pmt.SmoothLine):
    name, desc = 'smooth_line', 'SmoothLine'
    
class Traj10(pmt.SmoothBackAndForth):
    name, desc = 'smooth_back_and_forth', 'smooth back and forth'

class Traj11(pmt.CircleWithIntro):
    name, desc = 'circle_with_intro', 'circle with intro'

class Traj12(pmt.CircleWithIntro):
    name, desc = 'circle_with_intro_slow', 'circle with intro slow'
    def __init__(self): pmt.CircleWithIntro.__init__(self, c=[0, 0, -0.5], r=1., v=1., dt_intro=5., dt_stay=0.5)

class Traj13(pmt.CompositeTraj):
    name, desc = 'oval_with_intro', 'oval with intro'
    def __init__(self):
        oval = pmt.Oval(l=1, r=1, v=2.5, z=-0.25)
        Y0, Y1 = [-2.5, -1, 0, 0], oval.get(0.)
        takeoff = pmt.SmoothLine(Y0, Y1, duration=1.)
        Y2 =  [2.5, -1, 0, 0]
        landing = pmt.SmoothLine(oval.get(oval.duration), Y0, duration=4.)
        #goback = pmt.Line(Y2[:3], Y0[:3], v=1., psi=0)
        wait = pmt.Cst(Y0, duration=2.)
        steps = [takeoff, oval, landing, wait]
        pmt.CompositeTraj.__init__(self, steps)
        
trajectories = {T.name: (T.desc, T) for T in [Traj1, Traj2, Traj3, Traj4, Traj5, Traj6, Traj7, Traj8, Traj9, Traj10, Traj11, Traj12, Traj13]}

def get(traj_name):
    return trajectories[traj_name][1](), trajectories[traj_name][0]

def list():
    names = ['{}: {}'.format(k,v[0]) for k,v in sorted(trajectories.iteritems())]
    return names
