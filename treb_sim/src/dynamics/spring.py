'''springs'''

import numpy as np
from PyQt4 import QtGui, QtCore
from .misc import length_

DRAW_FORCE_SCALE = 400e3

class Spring(object):
    def __init__(self, sim=None, name=""):
        if(sim):
            self.sim = sim
            self.name = name
            self.idx = len(sim.springs)
            sim.springs.append(self)

class NailSpring(Spring):
    '''spring between object position and world position'''
    def __init__(self, sim, name, obj, xobj, x_world,
                 spring_constant, damping_constant=0):
        super().__init__(sim, name)
        self.frames = [obj.frame]
        self.frame = obj.frame
        self.x_frame = obj.obj2frame(np.mat(xobj).T)
        self.x_world = np.mat(x_world).T
        self.spring_constant = spring_constant
        self.damping_constant = damping_constant
        self.dim=2
        
    def displacement(self):
        '''return spring displacement [meters]'''
        return length_(self.frame.frame2world(self.x_frame) - self.x_world)
        

    def force(self):
        '''return column vector of [force_x, force_y, torque] at current time'''
        C    = self.frame.frame2world(self.x_frame) - self.x_world
        Cdot = self.frame.frame2worldv(self.x_frame)
        # torque is cross product (done in world space) of
        #  'arm' vector from frame's cg to spring point and
        #  force vector
        arm_vec = self.frame.frame2world(self.x_frame) - self.frame.origin
        force_vec = -(C * self.spring_constant +
                      Cdot * self.damping_constant).A1
        torque = np.asscalar(arm_vec[0]*force_vec[1] - arm_vec[1]*force_vec[0])
        return np.mat([force_vec[0], force_vec[1], torque]).T

    def draw(self, scene, draw_forces):
        '''draw the spring in a PyQt4 scene'''
        x_frame = self.frame.frame2world(self.x_frame).A1
        x_world = self.x_world
        pen = QtGui.QPen()
        pen.setColor(QtGui.QColor(128,128,128))
        pen.setWidth(0)

        scene.addLine(x_frame[0]-0.1, x_frame[1],
                      x_frame[0]+0.1, x_frame[1], pen)
        scene.addLine(x_frame[0], x_frame[1]-0.1,
                      x_frame[0], x_frame[1]+0.1, pen)
        scene.addEllipse(QtCore.QRectF(
                    QtCore.QPointF(x_world[0]-0.05, x_world[1]-0.05),
                    QtCore.QPointF(x_world[0]+0.05, x_world[1]+0.05)), pen)

        # draw force
        if (draw_forces):
            pen.setWidth(0.1)
            force = self.force()
            scene.addLine(x_frame[0],
                          x_frame[1],
                          x_frame[0]+force[0]/DRAW_FORCE_SCALE,
                          x_frame[1]+force[1]/DRAW_FORCE_SCALE, pen)

