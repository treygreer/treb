import numpy
from numpy import *
import pdb
from misc import length

DrawForceScale = 400e3

class Constraint:
    def __init__(self, sim=None, name=""):
        if(sim):
            self.sim = sim
            self.name = name
            self.idx = len(sim.constraints)
            sim.constraints.append(self)
        self.enabled = True

    def reset(self):
        pass

    def Fvec(self, i=0):
        forces = self.sim.F[:, self.idx, i]
        return numpy.sqrt(forces[:,0]**2 + forces[:,1]**2)

class Nail(Constraint):
    def __init__(self, sim, name, obj, xobj, xworld):
        Constraint.__init__(self, sim, name)
        self.frames = [obj.frame]
        self.frame = obj.frame
        self.xframe = obj.obj2frame(mat(xobj).T)
        self.xworld = mat(xworld).T
        self.dim=2

    def eval(self):
        C    = self.frame.frame2world(self.xframe) - self.xworld
        #if length(C) > 1e-3:
        #    print "nail error:  name=", self.name, "C=", C.T
        Cdot = self.frame.frame2worldv(self.xframe)
        # C = x + R*xframe - xworld
        # C0 = x0 + cos(theta)*xframe0 - sin(theta)*xframe1 - xworld0
        # C1 = x1 + sin(theta)*xframe0 + cos(theta)*xframe1 - xworld1
        j = hstack([eye(2), self.frame.dxdtheta(self.xframe)])
        # Cdot   = v + Rdot*xframe
        # Cdot0 = v0 + omega*(-sin(theta)*xframe0 - cos(theta)*xframe1)
        # Cdot1 = v1 + omega*( cos(theta)*xframe0 - sin(theta)*xframe1)
        jdot = hstack([zeros([2,2]), self.frame.dvdtheta(self.xframe)])
        return ({'C': C,
                 'Cdot': Cdot,
                 'blocks': [{'frame':  self.frame,
                             'j':    j,
                             'jdot': jdot}]})
    def draw(self, cr, drawForces):
        if (self.enabled):
            #print "Nail Force=", self.Force.A1
            xo = self.frame.frame2world(self.xframe).A1
            cr.set_source_rgb(0.5,0.5,0.5)
            pix,foo = cr.device_to_user_distance(1.0,1.0)
            cr.set_line_width(0.5*pix)
            cr.move_to(xo[0]-0.1, xo[1])
            cr.line_to(xo[0]+0.1, xo[1])
            cr.move_to(xo[0], xo[1]-0.1)
            cr.line_to(xo[0], xo[1]+0.1)
            cr.stroke()
            cr.arc(self.xworld.A1[0], self.xworld.A1[1], 0.05, 0, 2.0*pi)
            cr.stroke()
            # draw force
            if (drawForces):
                cr.set_line_width(4*pix)
                cr.move_to(xo[0],xo[1])
                f = [self.forces[0][0,0]/DrawForceScale,
                     self.forces[0][1,0]/DrawForceScale]
                cr.line_to(xo[0]+f[0],
                           xo[1]+f[1])
                cr.stroke()

class Pin(Constraint):
    def __init__(self, sim, name,
                 obj0, xobj0,
                 obj1, xobj1):
        Constraint.__init__(self, sim, name)
        self.frames = [obj0.frame, obj1.frame]
        self.frame0,self.frame1   = obj0.frame, obj1.frame
        self.xframe0 = obj0.obj2frame(mat(xobj0).T)
        self.xframe1 = obj1.obj2frame(mat(xobj1).T)
        self.dim=2
    def eval(self):
        C    = (self.frame0.frame2world(self.xframe0) -
                self.frame1.frame2world(self.xframe1))
        if length(C) > 1e-3:
            print "pin C=", C
            raise ValueError
        Cdot = (self.frame0.frame2worldv(self.xframe0) -
                self.frame1.frame2worldv(self.xframe1))
        j =    [ hstack([eye(2), self.frame0.dxdtheta(self.xframe0)]),
                -hstack([eye(2), self.frame1.dxdtheta(self.xframe1)])]
        jdot = [ hstack([zeros([2,2]), self.frame0.dvdtheta(self.xframe0)]),
                -hstack([zeros([2,2]), self.frame1.dvdtheta(self.xframe1)])]
        return ({'C': C,
                 'Cdot': Cdot,
                 'blocks': [ {'frame': self.frame0, 'j': j[0], 'jdot': jdot[0]},
                             {'frame': self.frame1, 'j': j[1], 'jdot': jdot[1]} ]})
    def draw(self, cr, drawForces):
        if (self.enabled):
            #print "Pin Force=", self.Force.A1
            xo0 = self.frame0.frame2world(self.xframe0).A1
            xo1 = self.frame1.frame2world(self.xframe1).A1
            cr.set_source_rgb(0.5,0.5,0.5)
            pix,foo = cr.device_to_user_distance(1.0,1.0)
            cr.set_line_width(0.5*pix)
            cr.move_to(xo0[0]-0.1, xo0[1])
            cr.line_to(xo0[0]+0.1, xo0[1])
            cr.move_to(xo0[0], xo0[1]-0.1)
            cr.line_to(xo0[0], xo0[1]+0.1)
            cr.stroke()
            cr.arc(xo1[0], xo1[1], 0.05, 0, 2.0*pi)
            cr.stroke()
            # draw force
            if (drawForces):
                cr.set_line_width(4*pix)
                cr.move_to(xo0[0],xo0[1])
                f = [self.forces[0][0,0]/DrawForceScale,
                     self.forces[0][1,0]/DrawForceScale]
                cr.line_to(xo0[0]+f[0],
                           xo0[1]+f[1])
                cr.stroke()

class Rod(Constraint):
    def __init__(self, sim, name,
                 obj0, xobj0, obj1, xobj1, length):
        Constraint.__init__(self, sim, name)
        self.frames = [obj0.frame, obj1.frame]
        self.frame0,self.frame1   = obj0.frame, obj1.frame
        self.xframe0 = obj0.obj2frame(mat(xobj0).T)
        self.xframe1 = obj1.obj2frame(mat(xobj1).T)
        self.length = length
        self.dim = 1
    def eval(self):
        xvec =  (self.frame1.frame2world(self.xframe1) -
                 self.frame0.frame2world(self.xframe0))
        vvec =  (self.frame1.frame2worldv(self.xframe1) -
                 self.frame0.frame2worldv(self.xframe0))
        C    = xvec.T*xvec - self.length*self.length
        if length(C) > 1e-3:
            print "rod error:  "
            print "  name=", self.name, "rodC=", C, "xvec=", xvec.T
            #raise ValueError
        Cdot = 2*xvec.T*vvec
        j    = [-2*hstack([xvec[0],
                           xvec[1],
                           self.frame0.dxdtheta(self.xframe0).T * xvec]),
                 2*hstack([xvec[0],
                           xvec[1],
                           self.frame1.dxdtheta(self.xframe1).T * xvec])]
        jdot = [-2*hstack([vvec[0],
                           vvec[1],
                           self.frame0.dvdtheta(self.xframe0).T * xvec +
                           self.frame0.dxdtheta(self.xframe0).T * vvec]),
                 2*hstack([vvec[0],
                           vvec[1],
                           self.frame1.dvdtheta(self.xframe1).T * xvec +
                           self.frame1.dxdtheta(self.xframe1).T * vvec])]
        return ({'C': C,
                 'Cdot': Cdot,
                 'blocks': [ {'frame': self.frame0, 'j': j[0], 'jdot': jdot[0]},
                             {'frame': self.frame1, 'j': j[1], 'jdot': jdot[1]} ]})
    def draw(self, cr, drawForces):
        if (self.enabled):
            #print "Rod Force=", self.Force.A1
            xo0 = self.frame0.frame2world(self.xframe0).A1
            xo1 = self.frame1.frame2world(self.xframe1).A1
            cr.set_source_rgb(0.5,0.5,0.5)
            pix,foo = cr.device_to_user_distance(1.0,1.0)
            cr.set_line_width(0.5*pix)
            cr.move_to(xo0[0], xo0[1])
            cr.line_to(xo1[0], xo1[1])
            cr.stroke()
            # draw force
            if (drawForces):
                cr.set_line_width(4*pix)
                cr.move_to(xo0[0],xo0[1])
                f = [self.forces[0][0,0]/DrawForceScale,
                     self.forces[0][1,0]/DrawForceScale]
                cr.line_to(xo0[0]+f[0],
                           xo0[1]+f[1])
                cr.stroke()
                cr.move_to(xo1[0],xo1[1])
                f = [self.forces[1][0,0]/DrawForceScale,
                     self.forces[1][1,0]/DrawForceScale]
                cr.line_to(xo1[0]+f[0],
                           xo1[1]+f[1])
                cr.stroke()

class Shelf(Nail):
    def __init__(self, sim, name, obj, xobj, height):
        Constraint.__init__(self, sim, name)
        self.frames = [obj.frame]
        self.frame = obj.frame
        self.xframe = obj.obj2frame(mat(xobj).T)
        self.xworld = mat([0, height]).T
        self.dim=1
    def eval(self):
        # evaluate as though a Nail, then pick out second (y) row of results
        results = Nail.eval(self)
        return ({'C': results['C'][1,:],
                 'Cdot': results['Cdot'][1,:],
                 'blocks': [ {'frame' : self.frame,
                              'j'   : results['blocks'][0]['j'][1,:],
                              'jdot': results['blocks'][0]['jdot'][1,:]}]})

    def draw(self, cr, drawForces):
        height = self.xworld.A1[1]
        xo = self.frame.frame2world(self.xframe).A1
        if (self.enabled):
            cr.set_source_rgb(0.0,0.0,0.0)
        else:
            cr.set_source_rgb(0.2,0.2,0.2)
        pix,foo = cr.device_to_user_distance(1.0,1.0)
        cr.set_line_width(1.0*pix)
        cr.move_to(-100.0, height)
        cr.line_to(100.0, height)
        cr.stroke()
#        cr.arc(xo[0], xo[1], 0.25, 0, 2.0*pi)
#        cr.stroke()

class Angle(Constraint):
    def __init__(self, sim, name, obj, theta):
        Constraint.__init__(self, sim, name)
        self.frames = [obj.frame]
        self.frame = obj.frame
        self.theta = obj.frame.theta + theta
        self.dim=1
    def eval(self):
        return ({'C':    mat([self.frame.theta - self.theta]),
                 'Cdot': mat([self.frame.omega]),
                 'blocks': [ {'frame' : self.frame,
                              'j'   : mat([0, 0, 1]),
                              'jdot': mat([0, 0, 0])} ]})
    def draw(self, cr, drawForces):
        pass

