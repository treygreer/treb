''' a bit more in the comment...
'''

import dynamics.simulation
from dynamics.frame import Frame
from dynamics.spring import NailSpring
from dynamics.object import Rectangle, Circle, Beam
from dynamics.constraint import Nail, Rod, Pin, Shelf
from dynamics.animation import Animation

from dynamics.constants import foot2meter, inch2meter, meter2foot
from dynamics.misc import length_, rot2radians, radians2rot
from dynamics.constants import lb2kgram, kgram2lb, newton2lb
from dynamics.constants import pine_density, steel_density

import scipy
import scipy.interpolate
import numpy as np
from math import pi, sin, cos, sqrt, acos, atan2
from scipy.optimize.minpack import fsolve

scipy.set_printoptions(precision=5, linewidth=200)

def continue_sim(sim, time, y):
    "continue simulation?"
    return True

sim_duration = 2.0
time_step = 0.001
debug = True

sim = dynamics.simulation.Simulation(max_time=sim_duration,
                                         time_step=time_step)
sim.debug=debug
hinge_pos = (0.0, 0.0)
cw_mass = lb2kgram(2000.)
ramp_length = foot2meter(10.)

sim.weightFrame=dynamics.frame.Frame(sim, "weight", theta=45*pi/180, origin=hinge_pos)
sim.weightFrame.ramp = dynamics.object.Rectangle(sim.weightFrame, l=ramp_length, w=inch2meter(4),
                                                 mass=0, color=(0.3,0.5,0.2),
                                                 origin = (ramp_length/2,0))
sim.weightFrame.cw = dynamics.object.Rectangle(sim.weightFrame, l=foot2meter(2.6), w=foot2meter(2.6),
                                               color=(0.3,0.5,0.2),
                                               mass=cw_mass,
                                               origin = (ramp_length,0))

# initialize frames
for frame in sim.frames:
    frame.init()

# define constraints
sim.hinge = Nail(sim, "hinge",
                 obj=sim.weightFrame.ramp,
                 xobj=(-ramp_length/2, 0.0),
                 xworld=(0,0))

print( "    running simulation")
sim.run(continue_sim, debug=debug)
print ("    done.")


def Y2range(sim, Y, with_air_friction=True):
    if (len(Y.shape)==1):
        Y = Y.reshape([1,len(Y)])
    idx = sim.weightFrame.idx
    x0 = Y[:,6*idx]
    y0 = Y[:,6*idx+1]
    vx0 = Y[:,6*idx+3]
    vy0 = Y[:,6*idx+4]
    
    tof = 2.0 * vy0 / scipy.constants.g
    tof[tof<0.0] = 0.0
    return (tof*vx0)

#anim=Animation(sim, Y2range)
