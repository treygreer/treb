from numpy import mat, array, zeros, copy, asarray
from scipy.integrate import ode
import pdb
from math import pi, sqrt
from dynamics.constants import air_density, inch2meter, lb2kgram
from scipy.constants import g

class Flight:
    def __init__(self,
                 mass=lb2kgram(8.0),
                 area=pi*inch2meter(4.0)**2,
                 Cd=0.2,
                 max_time=100, time_step=0.01,
                 debug=False):
        self.max_time  = max_time
        self.time_step = time_step
        self.Cd = Cd
        self.mass = mass
        self.area = area
        self.debug = debug
        self.num_steps = round(max_time / time_step) + 1

    def deriv(self, time, Y):
        x = Y[0]
        y = Y[1]
        vx = Y[2]
        vy = Y[3]

        # calculate air friction vector
        foo = air_density * self.Cd * self.area * sqrt(vx**2.0 + vy**2.0) / 2.0
        Fair_x = -vx * foo
        Fair_y = -vy * foo
        
        # calculate acceleration vector
        vdot_x = Fair_x / self.mass
        vdot_y = Fair_y / self.mass - g

        return(array([vx, vy, vdot_x, vdot_y]))

    def run(self, start, velocity):
        # state data structures
        self.t=zeros(self.num_steps)
        self.Y=zeros([self.num_steps, 4])
        self.Y[0] = array([start[0],start[1],velocity[0],velocity[1]])

        r = ode(self.deriv)
        r.set_integrator('vode', rtol=1e-9, atol=1e-6)
        r.set_initial_value(self.Y[0],t=0.0)

        time_idx=1
        while r.successful() and r.t < self.max_time:
            r.integrate(time_idx*self.time_step)
            self.t[time_idx] = r.t
            self.Y[time_idx] = r.y
            time_idx = time_idx + 1

            if self.debug:
                print ("time=%10.3f: x=%g y=%g vx=%g vy=%g" % \
                      (r.t, r.y[0], r.y[1], r.y[2], r.y[3]))

            if r.y[1]<0:
                break
            
        # simulation done: truncate arrays 
        self.Y.resize([time_idx, self.Y.shape[1]])
        self.t.resize([time_idx])

    def range(self):
        return self.Y[-1,0]
