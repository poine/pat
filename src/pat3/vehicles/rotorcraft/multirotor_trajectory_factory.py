import numpy as np

#
# A collection of multirotor trajectories
#

import pat3.trajectory_1D as p3_t1D
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
        r, v = 1., 2.; om = v/r; alpha0 = 0
        psit = p3_t1D.AffineOne(om, alpha0)
        pmt.Circle.__init__(self, [0, 0, -0.25], r, v, zt=None, psit=psit)
register(Traj3)


class Traj4(pmt.Circle):
    name, desc = 'circle3', 'circle r=2 v=4'
    def __init__(self): pmt.Circle.__init__(self, [0, 0, -0.5], r=1., v=4., psit=p3_t1D.AffineOne(1, 0))

class Traj5(pmt.Circle):
    name, desc = 'circle4', 'circle r=2 v=4'
    def __init__(self):
        r, v = 1., 4.; om = v/r
        pmt.Circle.__init__(self, c=[0, 0, -0.5], r=1., v=4., alpha0=0., dalpha=2*np.pi, psit=p3_t1D.SinOne(om=om))

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
# class Traj41(pmt.RefModTraj):
#     name, desc = 'refmod1', 'reference model test 1'
#     def __init__(self):
#         pmt.RefModTraj.__init__(self, [0, 0, -0.25], [1, 1, -0.25], v=0.1)
# register(Traj41)


##
# space indexed circle 1
#
import pat3.vehicles.rotorcraft.multirotor_trajectory_dev as trj_dev
class Traj42(trj_dev.SpaceIndexedTraj):
    name, desc = 'sic1', 'space index circle'
    def __init__(self, duration=10.):
        r, v = 1.5, 2.; om = v/r; alpha0 = 0
        psit = p3_t1D.CstOne(0.)
        #psit = pmt.AffineOne(om, alpha0+np.pi/2)
        #psit = pmt.PolynomialOne([-np.pi/2,0,0,0,0], [3*np.pi/2, 0, 0, 0, 0], duration=1.)
        zt = p3_t1D.CstOne(-1.5)
        straj = trj_dev.SpaceCircle(r=r, c=[-1.5,0], alpha0=0, dalpha=2*np.pi, ztraj=zt, psitraj=psit)
        dtraj = p3_t1D.PolynomialOne([0,0,0,0,0], [1, 0, 0, 0, 0], duration=duration)
        trj_dev.SpaceIndexedTraj.__init__(self,straj, dtraj)
register(Traj42)


##
#
#
class Traj43(trj_dev.SpaceIndexedTraj):
    name, desc = 'siwp1', 'space indexed waypoints'
    def __init__(self, duration=7., dtraj=None):
        if 1:
            straj = trj_dev.SpaceWaypoints([[0, 0, -1.5],
                                            [2, 1, -2.5],
                                            [0, 2, -1.2],
                                            [2, 3, -2.5],
                                            [0, 4, -1.5]]) # , [-2, 2, -2.]
        else:
            straj = trj_dev.SpaceWaypoints([[0, 0, -0.1],
                                            [2, 0, -0.6],
                                            [2, 2, -0.9],
                                            [0, 2, -0.6],
                                            [0, 0, -0.1],])
        #dtraj = pmt.AffineOne(1./duration,0., duration)
        dtraj = dtraj or p3_t1D.PolynomialOne([0,0,0,0,0], [1,0,0,0,0], duration)
        
        trj_dev.SpaceIndexedTraj.__init__(self,straj, dtraj)
register(Traj43)


##
#  quad4d rebooted
#
class Traj12_1(pmt.SmoothBackAndForth):
    name, desc = 'sbf', 'smooth back and forth'
    def __init__(self): pmt.SmoothBackAndForth.__init__(self, x0=[0, 0, -1.5, 0], x1=[1, 0, -2, np.pi/2])
register(Traj12_1)

class Traj12_2(pmt.CircleWithIntro):
    name, desc = 'cis', 'circle with intro slow'
    def __init__(self): pmt.CircleWithIntro.__init__(self, c=[0, 0, -2.5], r=1.5, v=1., dt_intro=5., dt_stay=0.5)
register(Traj12_2)

class Traj12_3(pmt.CircleWithIntro):
    name, desc = 'cis1', 'circle with intro slow'
    def __init__(self): pmt.CircleWithIntro.__init__(self, Y0=[0, 1, -1.5, 0], c=[0, 1, -2.5],
                                                     r=1.5, v=1., dt_intro=5., dt_stay=0.5, psit=p3_t1D.CstOne(0.))
register(Traj12_3)

class Traj12_4(pmt.CircleWithIntro):
    name, desc = 'cis2', 'circle with intro slow'
    def __init__(self): pmt.CircleWithIntro.__init__(self, Y0=[0, -1, -1.5, 0], c=[0, -1, -2.5],
                                                     r=-1.5, v=1., dt_intro=5., dt_stay=0.5, psit=p3_t1D.CstOne(0.))
register(Traj12_4)


#
# La sphere de gautier
#

class Sphere0:
    name, desc = 'sph0', 'sphere0'
    def __init__(self, c=[0, 0, -2.], r=1.5, v=2., psi=None):
        self.c, self.r, self.v = np.asarray(c), r, v # center, radius, velocity
        self.omega1, self.omega2 = 1, 0.1 #self.v/self.r
        self.t0, self.duration = 0, 80

        
    def reset(self, t0):
        self.t0 = t0

    def get(self, t):
        dt = t-self.t0
        alpha = self.omega1*(dt)# + self.alpha0
        beta = self.omega2*(dt)# + self.alpha0
        rca, rsa = np.abs(self.r)*np.cos(alpha), np.abs(self.r)*np.sin(alpha) 
        cb, sb = np.cos(beta), np.sin(beta)
        Yc = np.zeros((5,4))
        #Yc[0,:pmt._psi] = self.c[:pmt._psi] + [rca*cb, rsa, rca*sb]
        A = np.array([[cb, 0, -sb],[0, 1, 0],[sb, 0, cb]])
        B = np.array([rca, rsa, 0])
        Yc[0,:pmt._z+1] = self.c + A@B
        alpha_d, beta_d = self.omega1, self.omega2
        Ad = -beta_d*np.array([[sb, 0, cb],[0, 0, 0],[-cb, 0, sb]])
        Bd =  alpha_d * np.array([-rsa, rca, 0])
        Yc[1,:pmt._z+1] = Ad@B+A@Bd
        alpha_dd, beta_dd = 0, 0
        Add = -beta_dd*np.array([[sb, 0, cb],[0, 0, 0],[-cb, 0, sb]])-beta_d**2*np.array([[cb, 0, -sb],[0, 0, 0],[sb, 0, cb]])
        Bdd = alpha_dd*np.array([-rsa, rca, 0])-alpha_d**2*np.array([rca, rsa, 0])
        Yc[2,:pmt._z+1] = Add@B + 2*Ad@Bd + A@Bdd
        
        # pointing center
        #Yc[0,pmt._psi] = np.arctan2(Yc[0,pmt._y], Yc[0,pmt._x])
        om3 = 0.25
        #Yc[0,pmt._psi] = om3*dt
        #Yc[1,pmt._psi] = om3
        Yc[0,pmt._psi] =         np.sin(om3*dt)
        Yc[1,pmt._psi] =  om3   *np.cos(om3*dt)
        Yc[2,pmt._psi] = -om3**2*np.sin(om3*dt)
        Yc[3,pmt._psi] = -om3**3*np.cos(om3*dt)
        Yc[4,pmt._psi] =  om3**4*np.sin(om3*dt)
        
        return Yc.T
register(Sphere0)

class Donut0:
    name, desc = 'q4d_dnt1', 'quad4d rebooted: donut'
    def __init__(self, c=[0, 0, -2.], r=1., r2=1., v=4., psi=None, duration=80.):
        self.c, self.r, self.r2, self.v = np.asarray(c), r, r2, v # center, radius, velocity
        self.omega1, self.omega2 = 1, 0.1 #self.v/self.r
        self.t0, self.duration = 0, duration

    def reset(self, t0):
        self.t0 = t0

    def get(self, t):
        dt = t-self.t0
        alpha, beta = self.omega1*dt, self.omega2*dt
        rca, rsa = np.abs(self.r)*np.cos(alpha), np.abs(self.r)*np.sin(alpha)
        cb, sb = np.cos(beta), np.sin(beta)
        c1 = self.c + [-self.r2*sb, self.r2*cb, 0]
        A = np.array([[cb, -sb, 0],[sb, cb, 0],[0, 0, 1]])
        B = np.array([0, rsa, rca])
        Yc = np.zeros((5,4))
        Yc[0,:pmt._z+1] = c1 + A@B
        #cbd, sbd = self.omega2
        c1d = [-self.omega2*self.r2*cb, -self.omega2*self.r2*sb, 0]
        Ad = self.omega2 * np.array([[-sb, -cb, 0],[cb, -sb, 0],[0, 0, 1]])
        Bd = self.omega1 * np.array([0, rca, -rsa])
        Yc[1,:pmt._z+1] = c1d + Ad@B + A@Bd
        return Yc.T
register(Donut0)

class Donut1(pmt.CompositeTraj):
    name, desc = 'q4d_dnt2', 'quad4d rebooted: donut with intro'
    def __init__(self):
        Y0 = [0., 0, -1.5, 0.]
        d1 = Donut0(r=0.7, r2=1., duration=61.)
        Y1 = d1.get(0)#[:,0]
        Y2 = d1.get(d1.duration)
        steps = [pmt.SmoothLine(Y0, Y1, duration=2.),
                 d1,
                 pmt.SmoothLine(Y2, Y0, duration=2.),
                 pmt.Cst(Y0, duration=1.)]
        pmt.CompositeTraj.__init__(self, steps)


register(Donut1)


def print_available():
    print('available trajectories:')
    for n in list():
        print(' {}'.format(n))

def get(traj_name):
    return trajectories[traj_name][1](), trajectories[traj_name][0]

def list():
    names = ['{}: {}'.format(k,v[0]) for k,v in sorted(trajectories.items())]
    return names
