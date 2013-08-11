
import pdb
from numpy import mat,sqrt,array,sin,cos,hypot,hstack,linspace,trapz
from math import atan2,pi
import constants
from frame import *

class Object:
    "object class"
    def __init__(self, frame, cg=(0.0,0.0), origin=(0.0,0.0), theta=0.0):
        self.idx = len(frame.objects)
        frame.objects.append(self)
        self.frame = frame
        self.origin = mat(origin).T
        self.theta = theta
        self.update_working_vars()
        self.cg = mat(cg).T

    def update_working_vars(self):
        # r = [ cos(theta), sin(theta) ]
        self.r = radians2rot(self.theta)
        # R = |r0  -r1| = | cos  -sin |
        #     |r1   r0|   | sin   cos |
        self.R = rot2matrix(self.r)   # orientation matrix

    def draw(self, cr):
        cr.set_source_rgb(self.color[0], self.color[1], self.color[2])
        pix,foo = cr.device_to_user_distance(1.0,1.0)
        cr.set_line_width(0.5*pix)
        ### draw cg
        # xo = self.obj2world(self.cg)
        # cr.move_to(xo[0]-0.02, xo[1])
        # cr.line_to(xo[0]+0.02, xo[1])
        # cr.move_to(xo[0], xo[1]-0.02)
        # cr.line_to(xo[0], xo[1]+0.02)
        # cr.stroke()

    def Fext(self):
        "return world-space external force vector"
        # gravity only (no torque)
        return ([0, -scipy.constants.g*self.mass, 0])

    def PE(self):
        "return object potential energy"
        return (self.obj2world(self.cg)[1,0]*self.mass*scipy.constants.g)

    def KE(self):
        "return object kinetic energy"
        v = self.frame.frame2worldv(self.obj2frame(self.cg))
        return (0.5*self.mass*v.T*v +
                0.5*self.moment*self.frame.omega*self.frame.omega)

    def obj2frame(self, xobj):
        "transform from object space to frame space"
        return (self.origin + self.R*xobj)

    def frame2obj(self, xframe):
        "transform from frame space to object space"
        return (self.R.T*(xframe-self.origin))

    def obj2world(self, xobj):
        "transform from object space to world space"
        return (self.frame.frame2world(self.obj2frame(xobj)))

    def obj2worldv(self, xobj):
        "transform from object space position to world space velocity"
        return (self.frame.v + self.frame.Rdot*self.obj2frame(xobj))

    def dxdtheta(self, xobj):
        "given xobj, return d(xworld)/d(theta)"
        return (mat([[-self.frame.r[1], -self.frame.r[0]],
                     [ self.frame.r[0], -self.frame.r[1]]]) *
                self.obj2frame(xobj))

    def dvdtheta(self, xobj):
        "given xobj, return d(vworld)/d(theta)"
        return (self.frame.omega *
                mat([[-self.frame.r[0],  self.frame.r[1]],
                     [-self.frame.r[1], -self.frame.r[0]]]) *
                self.obj2frame(xobj))

    def world2obj(self, xworld):
        "transform from world space to object space"
        return (self.frame2obj(self.frame.world2frame(xworld)))

class Rectangle(Object):
    def __init__(self, frame, l, w, mass, color,
                 moment=None,
                 cg=(0.,0.), origin=(0.,0.), theta=0.0):
        Object.__init__(self, frame, cg, origin, theta)
        self.l,self.w,self.mass = l,w,mass
        self.color = color
        if moment is None:
            moment = (l*l + w*w) * mass/12.0
        self.moment = moment

        self.corners = mat([[-l/2.0, -w/2.0],
                            [ l/2.0, -w/2.0],
                            [ l/2.0,  w/2.0],
                            [-l/2.0,  w/2.0]]).T

    def draw(self, cr):
        Object.draw(self, cr)
        xw = self.obj2world(self.corners)
        pix,foo = cr.device_to_user_distance(1.0,1.0)
        cr.set_line_width(1.0*pix)
        cr.move_to(xw[0,0], xw[1,0])
        for i in [1,2,3]:
            cr.line_to(xw[0,i], xw[1,i])
        cr.close_path()
        cr.stroke()

class Circle(Object):
    def __init__(self, frame, radius, mass, color,
                 moment=None,
                 cg=(0.,0.), origin=(0.,0.), theta=0.0):
        Object.__init__(self, frame, cg, origin, theta)
        self.radius,self.mass = radius,mass
        self.color = color
        if moment is None:
            moment = mass * radius * radius / 2
        self.moment = moment

    def draw(self, cr):
        Object.draw(self, cr)
        pix,foo = cr.device_to_user_distance(1.0,1.0)
        cr.set_line_width(1.0*pix)
        center = self.obj2world(mat((0,0)).T)
        cr.arc(center[0,0], center[1,0], self.radius, 0, 2.0*pi)
        cr.stroke()

class Beam(Object):
    def __init__(self, frame,
                 x0, t0, d0,
                 x1, t1, d1,
                 density,
                 color,
                 origin=(0.,0.),
                 theta=0.0):
        Object.__init__(self, frame, origin=origin, theta=theta)
        self.x0,self.t0,self.d0 = x0,t0,d0
        self.x1,self.t1,self.d1 = x1,t1,d1
        self.color = color

        xvec = linspace(x0, x1, num=100)
        self.mass = trapz(self.t(xvec)*self.d(xvec)*density, xvec)
        xcg = trapz(xvec*self.t(xvec)*self.d(xvec)*density, xvec) / self.mass
        self.cg = mat((xcg,0.0)).T
        self.mass = abs(self.mass)
        self.moment = abs(trapz((xvec-xcg)**2.0 * self.t(xvec) * self.d(xvec) * density +
                                self.t(xvec) * self.d(xvec)**3.0 * density / 12.0,
                                xvec))

        self.corners = mat([[x0, -d0/2.0],
                            [x0,  d0/2.0],
                            [x1,  d1/2.0],
                            [x1, -d1/2.0]]).T

    def t(self, objx):
        "beam thickness (into screen) at objx"
        thick = self.t0*(self.x1-objx)/(self.x1-self.x0) + \
                self.t1*(objx-self.x0)/(self.x1-self.x0)
        if (any(thick <= 0.0)):
            print "thick = ", thick
            raise ValueError
        return thick

    def d(self, objx):
        "beam depth (height) at objx"
        depth = self.d0*(self.x1-objx)/(self.x1-self.x0) + \
                self.d1*(objx-self.x0)/(self.x1-self.x0)
        if (any(depth <= 0.0)):
            print "depth = ", depth
            raise ValueError
        return depth

    def draw(self, cr):
        Object.draw(self, cr)
        xw = self.obj2world(self.corners)
        pix,foo = cr.device_to_user_distance(1.0,1.0)
        cr.set_line_width(1.0*pix)
        cr.move_to(xw[0,0], xw[1,0])
        for i in [1,2,3]:
            cr.line_to(xw[0,i], xw[1,i])
        cr.close_path()
        cr.stroke()

