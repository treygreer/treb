import numpy as np
import pdb  # @UnusedImport
from .misc import length_
from .object import Object
from .frame import Frame
from math import pi
from PyQt4 import QtGui, QtCore

DrawForceScale = 400e3

class Constraint:
    def __init__(self, sim=None, name=""):
        if(sim):
            self.sim = sim
            self.name = name
            self.idx = len(sim.constraints)
            sim.constraints.append(self)
        self.enabled = True

    def Fvec(self, i=0):
        forces = self.sim.F[:, self.idx, i]
        return np.sqrt(forces[:,0]**2 + forces[:,1]**2)

class Nail(Constraint):
    def __init__(self, sim, name, obj, xobj, xworld):
        Constraint.__init__(self, sim, name)
        self.frames = [obj.frame]
        self.frame = obj.frame
        self.xframe = obj.obj2frame(np.mat(xobj).T)
        self.xworld = np.mat(xworld).T
        self.dim=2

    def eval(self):
        C    = self.frame.frame2world(self.xframe) - self.xworld
        if not isinstance(self, Shelf) and length_(C) > 1e-3:
            print("nail error:  name=", self.name, "C=", C.T)
            raise ValueError
        Cdot = self.frame.frame2worldv(self.xframe)
        j = np.hstack([np.eye(2), self.frame.dxdtheta(self.xframe)])
        jdot = np.hstack([np.zeros([2,2]), self.frame.dvdtheta(self.xframe)])
        return ({'C': C,
                 'Cdot': Cdot,
                 'blocks': [{'frame':  self.frame,
                             'j':    j,
                             'jdot': jdot}]})

    def draw(self, scene, drawForces):
        if self.enabled:
            xo = self.frame.frame2world(self.xframe).A1
            pen = QtGui.QPen()
            pen.setColor(QtGui.QColor(128,128,128))
            pen.setWidth(0)
        
            scene.addLine(xo[0]-0.1, xo[1],
                          xo[0]+0.1, xo[1], pen)
            scene.addLine(xo[0], xo[1]-0.1,
                          xo[0], xo[1]+0.1, pen)
            scene.addEllipse(QtCore.QRectF(
                        QtCore.QPointF(xo[0]-0.05, xo[1]-0.05),
                        QtCore.QPointF(xo[0]+0.05, xo[1]+0.05)), pen)

            # draw force
            if (drawForces):
                pen.setWidth(0.1)
                force = [self.forces[0].A1[0]/DrawForceScale,
                         self.forces[0].A1[1]/DrawForceScale]
                scene.addLine(xo[0],
                              xo[1],
                              xo[0]+force[0],
                              xo[1]+force[1], pen)

class Pin(Constraint):
    def __init__(self, sim, name,
                 obj0, xobj0,
                 obj1, xobj1=None):
        Constraint.__init__(self, sim, name)
        self.frame0 = obj0.frame
        self.xframe0 = obj0.obj2frame(np.mat(xobj0).T)
        if isinstance(obj1, Object):
            self.frame1 = obj1.frame
            self.xframe1 = obj1.obj2frame(np.mat(xobj1).T)
        else:
            assert isinstance(obj1, Frame)
            assert xobj1 is None
            self.frame1 = obj1
            self.xframe1 = \
                self.frame1.world2frame(self.frame0.frame2world(self.xframe0))
        self.frames = [self.frame0, self.frame1]
        self.dim=2
    def eval(self):
        C    = (self.frame0.frame2world(self.xframe0) -
                self.frame1.frame2world(self.xframe1))
        if length_(C) > 1e-3:
            print("pin {}: C = {}".format(self.name, C))
            raise ValueError
        Cdot = (self.frame0.frame2worldv(self.xframe0) -
                self.frame1.frame2worldv(self.xframe1))
        j =    [ np.hstack([np.eye(2), self.frame0.dxdtheta(self.xframe0)]),
                -np.hstack([np.eye(2), self.frame1.dxdtheta(self.xframe1)])]
        jdot = [ np.hstack([np.zeros([2,2]), self.frame0.dvdtheta(self.xframe0)]),
                -np.hstack([np.zeros([2,2]), self.frame1.dvdtheta(self.xframe1)])]
        return ({'C': C,
                 'Cdot': Cdot,
                 'blocks': [ {'frame': self.frame0, 'j': j[0], 'jdot': jdot[0]},
                             {'frame': self.frame1, 'j': j[1], 'jdot': jdot[1]} ]})

    def draw(self, scene, drawForces):
        if (self.enabled):
            pen = QtGui.QPen()
            pen.setColor(QtGui.QColor(128,128,128))
            pen.setWidth(0)
            #print "Pin Force=", self.Force.A1
            xo0 = self.frame0.frame2world(self.xframe0).A1
            xo1 = self.frame1.frame2world(self.xframe1).A1
            scene.addLine(xo0[0]-0.1, xo0[1],
                          xo0[0]+0.1, xo0[1], pen)
            scene.addLine(xo0[0], xo0[1]-0.1,
                          xo0[0], xo0[1]+0.1, pen)
            scene.addEllipse(QtCore.QRectF(
                        QtCore.QPointF(xo1[0]-0.05, xo1[1]-0.05),
                        QtCore.QPointF(xo1[0]+0.05, xo1[1]+0.05)), pen)

            # draw force
            if (drawForces):
                pen.setWidth(0.1)
                force = [self.forces[0].A1[0]/DrawForceScale,
                         self.forces[0].A1[1]/DrawForceScale]
                scene.addLine(xo0[0],
                              xo0[1],
                              xo0[0]+force[0],
                              xo0[1]+force[1], pen)

class Rod(Constraint):
    def __init__(self, sim, name,
                 obj0, xobj0, obj1, xobj1, length=None):
        Constraint.__init__(self, sim, name)
        self.frames = [obj0.frame, obj1.frame]
        self.frame0,self.frame1   = obj0.frame, obj1.frame
        self.xframe0 = obj0.obj2frame(np.mat(xobj0).T)
        self.xframe1 = obj1.obj2frame(np.mat(xobj1).T)
        self.length = length
        self.dim = 1
    def eval(self):
        xvec =  (self.frame1.frame2world(self.xframe1) -
                 self.frame0.frame2world(self.xframe0))
        vvec =  (self.frame1.frame2worldv(self.xframe1) -
                 self.frame0.frame2worldv(self.xframe0))
        if self.length is None:
            self.length = np.sqrt(np.asscalar(xvec.T * xvec))
        C    = xvec.T*xvec - self.length*self.length
        if length_(C) > 1e-3:
            print("rod error:  ")
            print("  name=", self.name, "rodC=", C, "xvec=", xvec.T)
            raise ValueError
        Cdot = 2*xvec.T*vvec
        j    = [-2*np.hstack([xvec[0],
                           xvec[1],
                           self.frame0.dxdtheta(self.xframe0).T * xvec]),
                 2*np.hstack([xvec[0],
                           xvec[1],
                           self.frame1.dxdtheta(self.xframe1).T * xvec])]
        jdot = [-2*np.hstack([vvec[0],
                           vvec[1],
                           self.frame0.dvdtheta(self.xframe0).T * xvec +
                           self.frame0.dxdtheta(self.xframe0).T * vvec]),
                 2*np.hstack([vvec[0],
                           vvec[1],
                           self.frame1.dvdtheta(self.xframe1).T * xvec +
                           self.frame1.dxdtheta(self.xframe1).T * vvec])]
        return ({'C': C,
                 'Cdot': Cdot,
                 'blocks': [ {'frame': self.frame0, 'j': j[0], 'jdot': jdot[0]},
                             {'frame': self.frame1, 'j': j[1], 'jdot': jdot[1]} ]})
    def draw(self, scene, drawForces):
        if (self.enabled):
            #print "Rod Force=", self.Force.A1
            pen = QtGui.QPen()
            pen.setColor(QtGui.QColor(128,128,128))
            pen.setWidth(0)
            xo0 = self.frame0.frame2world(self.xframe0).A1
            xo1 = self.frame1.frame2world(self.xframe1).A1
            scene.addLine(xo0[0], xo0[1],
                          xo1[0], xo1[1], pen)
            # draw force
            if (drawForces):
                pen.setWidth(0.1)
                force = [self.forces[0].A1[0]/DrawForceScale,
                         self.forces[0].A1[1]/DrawForceScale]
                scene.addLine(xo0[0],
                              xo0[1],
                              xo0[0]+force[0],
                              xo0[1]+force[1], pen)
                force = [self.forces[1].A1[0]/DrawForceScale,
                         self.forces[1].A1[1]/DrawForceScale]
                scene.addLine(xo1[0],
                              xo1[1],
                              xo1[0]+force[0],
                              xo1[1]+force[1], pen)

class Shelf(Nail):
    def __init__(self, sim, name, obj, xobj, height):
        Constraint.__init__(self, sim, name)
        self.frames = [obj.frame]
        self.frame = obj.frame
        self.xframe = obj.obj2frame(np.mat(xobj).T)
        self.xworld = np.mat([0, height]).T
        self.dim=1
    def eval(self):
        # evaluate as though a Nail, then pick out second (y) row of results
        results = Nail.eval(self)
        return ({'C': results['C'][1,:],
                 'Cdot': results['Cdot'][1,:],
                 'blocks': [ {'frame' : self.frame,
                              'j'   : results['blocks'][0]['j'][1,:],
                              'jdot': results['blocks'][0]['jdot'][1,:]}]})
    def draw(self, scene, drawForces):
        height = self.xworld.A1[1]
        pen = QtGui.QPen()
        pen.setWidth(0)
        if (self.enabled):
            pen.setColor(QtGui.QColor(0, 0, 0))
        else:
            pen.setColor(QtGui.QColor(51, 51, 51))
        scene.addLine(-100, height,
                      100, height, pen)

class Angle(Constraint):
    def __init__(self, sim, name, obj, theta):
        Constraint.__init__(self, sim, name)
        self.frames = [obj.frame]
        self.frame = obj.frame
        self.theta = obj.frame.theta + theta
        self.dim=1
    def eval(self):
        return ({'C':    np.mat([self.frame.theta - self.theta]),
                 'Cdot': np.mat([self.frame.omega]),
                 'blocks': [ {'frame' : self.frame,
                              'j'   : np.mat([0, 0, 1]),
                              'jdot': np.mat([0, 0, 0])} ]})
    def draw(self, scene, drawForces):
        pass

