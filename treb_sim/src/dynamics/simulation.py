import numpy as np
from numpy import mat, zeros
from scipy.integrate import ode

class Simulation:
    def __init__(self, max_time, time_step):
        self.max_time  = max_time
        self.time_step = time_step
        self.num_steps = round(max_time / time_step) + 1

        # numeric drift correction constants
        self.ks = 0 #1  # proportional restoring force
        self.kd = 0 #1  # differential restoring force

        self.frames = []
        self.constraints = []
        self.springs = []

        self.num_calls = 0

    def _Cdim(self):
        """Return the sum of the simulation's enabled constraint dimensions."""
        Cdim = 0
        for constraint in self.constraints:
            if constraint.enabled:
                Cdim = Cdim + constraint.dim
        return Cdim

    def deriv(self, time, y):
        """Given a time and state vector y, returns the derivative dy."""
        self.num_calls = self.num_calls+1

        # distribute state across frames for use in constraints
        self._scatter_state(y)

        # build velocity vector v
        v = mat(zeros([3*len(self.frames),1]))
        i=0
        for frame in self.frames:
            v[3*i:3*i+3,0] = mat(y[6*i+3:6*i+6]).T
            i=i+1

        # build constraint Jacobian and its time derivative
        Cdim = self._Cdim()
        J    = mat(zeros(shape=[Cdim, 3*len(self.frames)]))
        Jdot = mat(zeros(shape=[Cdim, 3*len(self.frames)]))
        C    = mat(zeros(shape=[Cdim,1]))
        Cdot = mat(zeros(shape=[Cdim,1]))
        row = 0
        for constraint in self.constraints:
            if constraint.enabled:
                ev = constraint.eval()
                C   [row:row+constraint.dim,0] = ev['C']
                Cdot[row:row+constraint.dim,0] = ev['Cdot']
                for block in ev['blocks']:
                    fidx = block['frame'].idx   # frame index
                    J   [row:row+constraint.dim,3*fidx:3*fidx+3] = block['j']
                    Jdot[row:row+constraint.dim,3*fidx:3*fidx+3] = block['jdot']
                row = row + constraint.dim

        Fspring = mat(zeros(shape=[3*len(self.frames), 1]))
        for spring in self.springs:
            idx = spring.frame.idx 
            Fspring[3*idx:3*(idx+1), 0] += spring.force()

        # solve linear constraint system
        A = J * self.Minv * J.T
        b = -(Jdot*v + J*self.Minv*(self.Fext+Fspring) + self.ks*C + self.kd*Cdot)
        lambda_ = A.I * b

        # write forces back into constraints and springs for 
        #   diagnostics and plotting
        row = 0   # row index of J
        for constraint in self.constraints:
            if constraint.enabled:
                constraint.lambda_ = lambda_[row:row+constraint.dim]
                Force = J.T[:,row:row+constraint.dim] * constraint.lambda_
                constraint.forces = []
                for frame in constraint.frames:
                    F = Force[3*frame.idx:3*frame.idx+3]
                    force = F[0:2,:]
                    constraint.forces.append(force)
                row = row + constraint.dim
            else:
                constraint.forces = None

        # calculate system accelerations
        F = J.T * lambda_
        vdot = self.Minv * (F + self.Fext + Fspring)

        # gather dy from frames and return
        dy   = zeros(6*len(self.frames))
        i=0
        for frame in self.frames:
            dy[6*i+0:6*i+3] = v.A   [3*i:3*i+3,0] 
            dy[6*i+3:6*i+6] = vdot.A[3*i:3*i+3,0]
            i=i+1

        return(dy)

    def _scatter_state(self, y):
        for frame in self.frames:
            frame.scatter_state(y)

    def _gather_state(self, y):
        for frame in self.frames:
            frame.gather_state(y)

    def run(self, continue_fun=None, debug=True):
        # build mass matrix M and external force vector Fext
        self.Minv = mat(zeros(shape=[3*len(self.frames), 3*len(self.frames)]))
        self.Fext = mat(zeros(shape=[3*len(self.frames), 1]))
        f=0
        for frame in self.frames:
            self.Minv[3*f:3*(f+1), 3*f:3*(f+1)] = frame.Minv()
            self.Fext[3*f:3*(f+1),0] = mat(frame.Fext()).T
            f=f+1

        # state data structure
        self.Y=zeros([self.num_steps, 6*len(self.frames)])

        # other saved data
        self.t=zeros(self.num_steps)
        self.F=zeros([self.num_steps, len(self.constraints), 2, 2])
        self.constraints_enabled=zeros([self.num_steps, len(self.constraints)])

        # gather initial state from frames
        self._gather_state(self.Y[0])
        c=0
        for constraint in self.constraints:
            self.constraints_enabled[0,c] = constraint.enabled
            c=c+1

        r = ode(self.deriv)
        r.set_integrator('vode', rtol=1e-9, atol=1e-6)
        r.set_initial_value(self.Y[0],t=0.0)

        time_idx=1
        while r.successful() and r.t < self.max_time:
            r.integrate(time_idx*self.time_step)
            self.t[time_idx] = r.t
            self.Y[time_idx] = r.y

            # save constraint forces
            c=0
            for constraint in self.constraints:
                self.constraints_enabled[time_idx,c] = constraint.enabled
                if constraint.enabled:
                    for f in range(len(constraint.frames)):
                        self.F[time_idx, c, f, :] = constraint.forces[f].T
                c=c+1

            time_idx = time_idx + 1

            if continue_fun and not continue_fun(self, r.t, r.y):
                break

        # simulation done: truncate arrays
        self.Y.resize([time_idx, self.Y.shape[1]])
        self.F.resize([time_idx] + list(self.F.shape[1:4]))
        self.t.resize([time_idx])
        self.constraints_enabled.resize([time_idx, self.constraints_enabled.shape[1]])
        
        # update minimum potential energies
        for frame in self.frames:
            frame.PEmin = np.min(frame.PEvec())
