import numpy as np

import pat3.vehicles.rotorcraft.multirotor_trajectory as pmt

trajectories = {}
def register(T): trajectories[T.name] = (T.desc, T)
##
#
#
class Traj1(pmt.Line):
    name, desc = 'l1', 'line v=0.1'
    def __init__(self): pmt.Line.__init__(self, [0, 0, -0.25], [1, 1, -0.25], v=0.1)
register(Traj1)

##
#
#
class Traj2(pmt.Circle):
    name, desc = 'c1', 'circle r=1 v=1, looking forward'
    def __init__(self):
        r, v = 1., 1.; om = v/r; alpha0 = 0
        #psit = pmt.AffineOne(om, alpha0+np.sign(r)*np.pi/2)
        #psit=pmt.SinOne(om=om)
        psit = None
        pmt.Circle.__init__(self, [0, 0, -0.25], r, v, zt=None, psit=psit)
register(Traj2)

##
#
#
class Traj3(pmt.Circle):
    name, desc = 'c2', 'circle r=1 v=4, looking at center'
    def __init__(self):
        r, v = 1., 4.; om = v/r; alpha0 = 0
        psit = pmt.AffineOne(om, alpha0)
        pmt.Circle.__init__(self, [0, 0, -0.25], r, v, zt=None, psit=psit)
register(Traj3)


class Traj4(pmt.Circle):
    name, desc = 'circle3', 'circle r=2 v=4'
    def __init__(self): pmt.Circle.__init__(self, [0, 0, -0.5], r=1., v=4., psit=pmt.AffineOne(1, 0))

class Traj5(pmt.Circle):
    name, desc = 'circle4', 'circle r=2 v=4'
    def __init__(self):
        r, v = 1., 4.; om = v/r
        pmt.Circle.__init__(self, c=[0, 0, -0.5], r=1., v=4., alpha0=0., dalpha=2*np.pi, psit=pmt.SinOne(om=om))

##
#  Oval
#
class Traj6(pmt.Oval):
    name, desc = 'oval1', 'oval'
    def __init__(self): pmt.Oval.__init__(self, l=1, r=1, v=3., z=-0.5)
register(Traj6)

class Traj7(pmt.DoubleOval):
    name, desc = 'doval1', 'double oval'
    def __init__(self): pmt.DoubleOval.__init__(self, l=2, r=1, v=4)
    
##
#  figure of eight
#
class Traj8(pmt.FigureOfEight):
    name, desc = 'foe1', 'figure of eight'
register(Traj8)


class Traj9(pmt.SmoothLine):
    name, desc = 'smooth_line', 'SmoothLine'
register(Traj9)
    
class Traj10(pmt.SmoothBackAndForth):
    name, desc = 'smooth_back_and_forth', 'smooth back and forth'


class Traj14(pmt.CompositeTraj):
    name, desc = 'smooth_back_and_forth2', 'SmoothBackAndForth2'
    def __init__(self):
        Y0, Y1 = [0., 0, -0.25, np.pi/2],  [0.5, 0, -0.25, np.pi/2]
        steps = [pmt.SmoothLine(Y0, Y1, duration=2.),
                 pmt.Cst(Y1, duration=1.),
                 pmt.SmoothLine(Y1, Y0, duration=2.),
                 pmt.Cst(Y0, duration=1.)]
        pmt.CompositeTraj.__init__(self, steps)

class Traj15(pmt.CompositeTraj):
    name, desc = 'smooth_back_and_forth3', 'SmoothBackAndForth3'
    def __init__(self):
        Y0, Y1 = [0., 0, -0.25, np.pi/2-np.deg2rad(30)],  [0., 0, -0.25, np.pi/2+np.deg2rad(30)]
        steps = [pmt.SmoothLine(Y0, Y1, duration=2.),
                 pmt.Cst(Y1, duration=1.),
                 pmt.SmoothLine(Y1, Y0, duration=2.),
                 pmt.Cst(Y0, duration=1.)]
        pmt.CompositeTraj.__init__(self, steps)

class Traj16(pmt.CompositeTraj):
    name, desc = 'smooth_back_and_forth4', 'SmoothBackAndForth4'
    def __init__(self):
        Y0, Y1 = [0., 0, -0.25, np.pi],  [0., 0, -0.5, np.pi]
        steps = [pmt.SmoothLine(Y0, Y1, duration=1.),
                 pmt.Cst(Y1, duration=1.),
                 pmt.SmoothLine(Y1, Y0, duration=1.),
                 pmt.Cst(Y0, duration=1.)]
        pmt.CompositeTraj.__init__(self, steps)        
##
#  figure of eight
#
class Traj17(pmt.CompositeTraj):
    name, desc = 'smooth_back_and_forth5', 'SmoothBackAndForth5'
    def __init__(self):
        Y0, Y1, Y2 = [0.5, 0.5, -0.25, np.pi],  [-0.5, 0.5, -0.25, 0],  [0., 0, -0.5, np.pi/2]
        steps = [pmt.SmoothLine(Y0, Y1, duration=1.),
                 pmt.Cst(Y1, duration=0.2),
                 pmt.SmoothLine(Y1, Y2, duration=1.),
                 pmt.Cst(Y2, duration=0.2),
                 pmt.SmoothLine(Y2, Y0, duration=1.),
                 pmt.Cst(Y0, duration=0.2)]
        pmt.CompositeTraj.__init__(self, steps)           
register(Traj17)


class Traj11(pmt.CircleWithIntro):
    name, desc = 'circle_with_intro', 'circle with intro'

class Traj12(pmt.CircleWithIntro):
    name, desc = 'circle_with_intro_slow', 'circle with intro slow'
    def __init__(self): pmt.CircleWithIntro.__init__(self, c=[0, 0, -0.5], r=1., v=1., dt_intro=5., dt_stay=0.5)

##
#  oval_with_intro
#
class Traj13(pmt.CompositeTraj):
    name, desc = 'oval_with_intro', 'oval with intro'
    def __init__(self):
        oval = pmt.Oval(l=1, r=1, v=2.5, z=-0.25)
        Y0, Y1 = [-2.5, -1, 0, np.pi/2], oval.get(0.)
        takeoff = pmt.SmoothLine(Y0, Y1, duration=1.)
        Y2 =  [2.5, -1, 0, 0]
        landing = pmt.SmoothLine(oval.get(oval.duration), Y0, duration=5.)
        #goback = pmt.Line(Y2[:3], Y0[:3], v=1., psi=0)
        wait = pmt.Cst(Y0, duration=2.)
        steps = [takeoff, oval, landing, wait]
        pmt.CompositeTraj.__init__(self, steps)
register(Traj13)

##
#  Smooth oval
#
class Traj18(pmt.CompositeTraj):
    name, desc = 'smooth_oval', 'smooth oval'
    def __init__(self, l=1., r=1, z=-0.25, v=3):
        c1, c2 = np.array([-l, 0, z]), np.array([l, 0, z])
        p1, p2 = np.array([-l, -r, z]), np.array([l, -r, z])
        eps_l, eps_alpha = 0.5, np.deg2rad(10)
        p1bis, p2bis = np.array([-l+eps_l, -r, z]), np.array([l-eps_l, -r, z])
        p3, p4 = np.array([l, r, z]), np.array([-l, r, z])
        circle1 = pmt.Circle(c2, r, v, -np.pi/2, np.pi)
        line1 = pmt.Line(p1bis, p2bis, v)
        circle2 = pmt.Circle(c1, r, v, np.pi/2, np.pi)
        Y0 = circle2.get(circle2.duration)
        Y1 = circle1.get(0)
        Y2 = circle1.get(circle1.duration)
        Y3 = circle2.get(0)
        steps = [#pmt.SmoothLine(circle2.get(circle2.duration), line1.get(0), duration=0.16),
                 #line1,
                 #pmt.SmoothLine(line1.get(line1.duration), Y1, duration=0.16),
                 pmt.SmoothLine(Y0, Y1, duration=0.7),
                 circle1,
                 pmt.SmoothLine(Y2, Y3, duration=0.7),
                 circle2]
        pmt.CompositeTraj.__init__(self, steps)
register(Traj18)


class SmoothOval2(pmt.CompositeTraj):
    name, desc = 'smooth_oval2', 'smooth oval2'
    def __init__(self, l=1.1, r=1, z=-0.25, v=3, dt_trans=0.203, dalpha=np.deg2rad(16), eps_l=0.33):
        c1, c2 = np.array([-l, 0, z]), np.array([l, 0, z])
        circle1 = pmt.Circle(c2, r, v, -np.pi/2+dalpha, np.pi-2*dalpha)
        circle2 = pmt.Circle(c1, r, v, np.pi/2+dalpha, np.pi-2*dalpha)
        l1s, l1e =  np.array([-l+eps_l, -r, z]), np.array([l-eps_l, -r, z])
        line1 = pmt.Line(l1s, l1e, v)
        l2s, l2e =  np.array([l-eps_l, r, z]), np.array([-l+eps_l, r, z])
        line2 = pmt.Line(l2s, l2e, v)
        trans1 = pmt.SmoothLine(line1.get(line1.duration), circle1.get(0), duration=dt_trans)
        trans2 = pmt.SmoothLine(circle1.get(circle1.duration), line2.get(0), duration=dt_trans)
        trans3 = pmt.SmoothLine(line2.get(line2.duration), circle2.get(0), duration=dt_trans)
        trans4 = pmt.SmoothLine(circle2.get(circle2.duration), line1.get(0), duration=dt_trans)
        steps = [line1, trans1, circle1, trans2, line2, trans3, circle2, trans4]
        pmt.CompositeTraj.__init__(self, steps)
        
register(SmoothOval2)     


##
#
#
class Traj19(pmt.Circle):
    name, desc = 'c42', 'circle r=30 v=15, looking forward'
    def __init__(self):
        r, v = 20., 12.; om = v/r; alpha0 = 0
        #psit = pmt.AffineOne(om, alpha0+np.sign(r)*np.pi/2)
        #psit=pmt.SinOne(om=om)
        psit = None
        pmt.Circle.__init__(self, [0, 0, -0.25], r, v, zt=None, psit=psit)
register(Traj19)



##
#
#
class Traj41(pmt.RefModTraj):
    name, desc = 'refmod1', 'reference model test 1'
    def __init__(self):
        pmt.RefModTraj.__init__(self, [0, 0, -0.25], [1, 1, -0.25], v=0.1)
register(Traj41)


##
#
#
import pat3.vehicles.rotorcraft.multirotor_trajectory_dev as trj_dev
class Traj42(trj_dev.SpaceIndexedTraj):
    name, desc = 'si1', 'space index circle'
    def __init__(self, duration=10.):
        r, v = 1.5, 2.; om = v/r; alpha0 = 0
        psit = pmt.CstOne(0.)
        #psit = pmt.AffineOne(om, alpha0+np.pi/2)
        #psit = pmt.PolynomialOne([-np.pi/2,0,0,0,0], [3*np.pi/2, 0, 0, 0, 0], duration=1.)
        straj = trj_dev.SpaceCircle(r=r, c=[2.,2.], alpha0=0, dalpha=2*np.pi, psitraj=psit)
        dtraj = pmt.PolynomialOne([0,0,0,0,0], [1, 0, 0, 0, 0], duration=duration)
        trj_dev.SpaceIndexedTraj.__init__(self,straj, dtraj)
register(Traj42)


##
#
#
class Traj43(trj_dev.SpaceIndexedTraj):
    name, desc = 'si2', 'space indexed waypoints'
    def __init__(self, duration=7., dtraj=None):
        if 1:
            straj = trj_dev.SpaceWaypoints([[0, 0, -0.1],
                                            [2, 1, -0.8],
                                            [0, 2, -0.4],
                                            [2, 3, -0.8],
                                            [0, 4, -0.1]])
        else:
            straj = trj_dev.SpaceWaypoints([[0, 0, -0.1],
                                            [2, 0, -0.6],
                                            [2, 2, -0.9],
                                            [0, 2, -0.6],
                                            [0, 0, -0.1],])
        #dtraj = pmt.AffineOne(1./duration,0., duration)
        dtraj = dtraj or pmt.PolynomialOne([0,0,0,0,0], [1,0,0,0,0], duration)
        
        trj_dev.SpaceIndexedTraj.__init__(self,straj, dtraj)
register(Traj43)


#trajectories = {T.name: (T.desc, T) for T in [Traj1, Traj2, Traj3, Traj4, Traj5, Traj6, Traj7, Traj8, Traj9, Traj10, Traj11, Traj12, Traj13, Traj14, Traj15, Traj16, Traj17, Traj18]}

def print_available():
    print('available trajectories:')
    for n in list():
        print(' {}'.format(n))

def get(traj_name):
    return trajectories[traj_name][1](), trajectories[traj_name][0]

def list():
    names = ['{}: {}'.format(k,v[0]) for k,v in sorted(trajectories.items())]
    return names
