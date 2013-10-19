import numpy as np
import pdb  # @UnusedImport
from .misc import length_
from math import pi
from PyQt4 import QtGui, QtCore

DrawForceScale = 400e3

class Spring:
    def __init__(self, sim=None, name=""):
        if(sim):
            self.sim = sim
            self.name = name
            self.idx = len(sim.springs)
            sim.springs.append(self)

    def Fvec(self, i=0):
        forces = self.sim.F[:, self.idx, i]
        return np.sqrt(forces[:,0]**2 + forces[:,1]**2)

class NailSpring(Spring):
    '''spring between object position and world position'''
    def __init__(self, sim, name, obj, xobj, xworld, 
                 spring_constant):
        super().__init__(sim, name)
        self.frames = [obj.frame]
        self.frame = obj.frame
        self.xframe = obj.obj2frame(np.mat(xobj).T)
        self.xworld = np.mat(xworld).T
        self.spring_constant = spring_constant
        self.dim=2

    def force(self):
        C    = self.frame.frame2world(self.xframe) - self.xworld
        # torque is cross product (done in world space) of 
        #  'arm' vector from frame's cg to spring point and 
        #  force vector 
        arm_vec = self.frame.frame2world(self.xframe) - self.frame.origin
        force_vec = (-C * self.spring_constant).A1
        torque = np.asscalar(arm_vec[0]*force_vec[1] - arm_vec[1]*force_vec[0])
        return np.mat([force_vec[0], force_vec[1], torque]).T

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

