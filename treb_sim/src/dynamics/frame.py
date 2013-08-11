import pdb
import numpy
import scipy
from scipy import constants
from scipy import mat,sqrt,array,sin,cos,hypot,hstack,zeros
from math import atan2,pi
from misc import *
import constants

class Frame:
    "frame class"
    def __init__(self, sim, name, origin=(0.0,0.0), theta=0.0):
        self.sim = sim
        self.name = name
        self.idx = len(sim.frames)
        sim.frames.append(self)
        self.origin = mat(origin).T
        self.theta = theta
        self.omega = 0.0
        self.v = mat((0.0,0.0)).T
        self.update_working_vars()
        self.objects = []
        self.PEmin = float('inf')

    def init(self, cg=None, mass=None, moment=None):
        if (cg or mass or moment):
            assert (cg and mass and moment)
            self.mass = mass
            self.cg = cg
        else:
            # compute frame's mass and center of gravity
            self.cg = mat((0.0,0.0)).T
            self.mass = 0.0
            for obj in self.objects:
                self.mass += obj.mass
                self.cg += obj.mass * obj.obj2frame(obj.cg)
            self.cg /= self.mass

        # The frame's center of gravity must be at its origin for the motion equations
        #  to work.  Offset all objects and the frame's initial position.
        self.origin = self.frame2world(self.cg)  # update frame origin in worldspace
        for obj in self.objects:
            obj.origin -= self.cg  # update object origin in framespace

        # compute the frame's moment based on new object positions about frame's cg
        # (frame origin)
        if moment:
            self.moment = moment
        else:
            self.moment = 0.0
            for obj in self.objects:
                self.moment += obj.moment + \
                               obj.mass * length(obj.obj2frame(obj.cg))**2

    def scatter_state(self, y):
        my_y = y[6*self.idx:6*self.idx+6]
        self.origin = mat(my_y[0:2]).T
        self.theta = my_y[2]
        self.v = mat(my_y[3:5]).T
        self.omega = my_y[5]
        self.update_working_vars()

    def update_working_vars(self):
        # r = [ cos(theta), sin(theta) ]
        self.r = radians2rot(self.theta)
        # R = |r0  -r1| = | cos  -sin |
        #     |r1   r0|   | sin   cos |
        self.R = rot2matrix(self.r)   # orientation matrix
        # Rdot = |0  -w|  x  |r0  -r1|  =  w*|-r1  -r0| = w*| -sin  -cos |
        #        |w   0|     |r1   r0|       | r0  -r1|     |  cos  -sin |
        self.Rdot = dual(self.omega) * self.R 

    def gather_state(self, y):
        y[6*self.idx:6*self.idx+6] = hstack((self.origin.A1,
                                             self.theta,
                                             self.v.A1,
                                             self.omega))

    def draw(self, cr):
        for obj in self.objects:
            obj.draw(cr)
        # draw cg
        xo = self.origin
        cr.set_source_rgb(0.0,0.0,0.0)
        pix,foo = cr.device_to_user_distance(1.0,1.0)
        cr.set_line_width(0.5*pix)
        cr.move_to(xo[0]-0.02, xo[1])
        cr.line_to(xo[0]+0.02, xo[1])
        cr.move_to(xo[0], xo[1]-0.02)
        cr.line_to(xo[0], xo[1]+0.02)
        cr.stroke()
            
    def Fext(self):
        "return external force vector"
        # gravity only (no torque)
        return ([0, -scipy.constants.g*self.mass, 0])

    def Minv(self):
        "return mass inverse matrix"
        return ([[1.0/self.mass, 0, 0],
                 [0, 1.0/self.mass, 0],
                 [0, 0, 1.0/self.moment]])
                 
    def Mdiag(self):
        "return mass matrix diagonal"
        return ([self.mass, self.mass, self.moment])

    def PE(self):
        "return potential energy at current time"
        return (self.origin[1,0]*self.mass*scipy.constants.g)

    def PEvec(self):
        "return potential energy vector"
        Y = self.sim.Y
        PE=zeros(Y.shape[0])
        for time_idx in range(Y.shape[0]):
            self.scatter_state(Y[time_idx, :])
            PE[time_idx] = self.origin[1,0]*self.mass*scipy.constants.g
        return(PE)

    def KE(self):
        "return kinetic energy at current time"
        return (0.5*self.mass*self.v.T*self.v +
                0.5*self.moment*self.omega*self.omega)

    def KEvec(self):
        "return kinetic energy vector"
        Y = self.sim.Y
        KE=zeros(Y.shape[0])
        for time_idx in range(Y.shape[0]):
            self.scatter_state(Y[time_idx, :])
            KE[time_idx] = (0.5*self.mass*self.v.T*self.v +
                            0.5*self.moment*self.omega*self.omega)
        return(KE)

    def frame2world(self, xframe):
        "transform from frame space to world space"
        return (self.origin + self.R*xframe)

    def world2frame(self, xworld):
        "transform from world space to frame space"
        return (self.R.T*(xworld-self.origin))

    def frame2worldv(self, xframe):
        "transform from frame space position to world space velocity"
        return (self.v + self.Rdot*xframe)

    def dxdtheta(self, xframe):
        "given xframe, return d(xworld)/d(theta)"
        return (mat([[-self.r[1], -self.r[0]],
                     [ self.r[0], -self.r[1]]]) * xframe)

    def dvdtheta(self, xframe):
        "given xframe, return d(vworld)/d(theta)"
        return (self.omega * mat([[-self.r[0],  self.r[1]],
                                  [-self.r[1], -self.r[0]]]) * xframe)

