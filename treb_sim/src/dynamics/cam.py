from numpy import mat,linspace,zeros,hstack,vstack,isscalar,diff,pi
from numpy import sqrt,sin,cos,cumsum,mod,sign,arctan2
from constraint import Constraint,Rod
from scipy import optimize
from scipy.interpolate import splprep, splev
from scipy.interpolate.fitpack2 import UnivariateSpline
import convex_hull
from object import length
import pdb

class Cam(Rod):
    def __init__(self, sim=None, name="", boundary=None,
                 cam_object=None, cam_center=(0,0), 
                 wt_object=None, wt_attach=(0,0),
                 length=None):
        Constraint.__init__(self, sim, name=name)
        self.frames = [cam_object.frame, wt_object.frame]
        self.cam_frame = cam_object.frame
        self.center = cam_object.obj2frame(mat(cam_center).T)
        self.camu = 0.0  # cam parameter in range [0,1]
        self.frame0 = self.cam_frame
        self.frame1, self.xframe1 = wt_object.frame, wt_object.obj2frame(mat(wt_attach).T)
        self.total_length = length
        self.dim = 1
        boundary = cam_object.obj2frame(boundary.T).T + self.center.T
        boundary = convex_hull.hull(boundary)
        # periodic cubic spline fit to boundary
        self.camtck,foo = splprep([boundary[:,0],boundary[:,1]],
                                  per=True, quiet=0, k=3, s=0)

        # build interpolation spline for line integral
        u=linspace(0,1,num=1000)
        v=splev(u, self.camtck)
        l = zeros(len(u))
        diffx = diff(v[0])
        diffy = diff(v[1])
        l = sqrt(diffx**2.0 + diffy**2.0)
        l = cumsum(hstack([0, l]))
        self.lspline = UnivariateSpline(u, l, s=0)
        self.camu = 0.5

    def tangent_point(self, u):
        "return cam's frame-space tangent point given cam parameter u"
        v=splev(mod(u, 1.0), self.camtck)
        return mat(vstack([v[0], v[1]]))

    def tangent_u(self, pivot, u0=0.0, disp=False):
        """Returns cam's parameter 'u' for the tangent condition given a
        world-space point 'pivot' on the tangent line.  Optional parameter
        'u0' is a starting guess for the cam parameter."""
        xpivot = self.cam_frame.world2frame(pivot)  # convert to cam's frame space

        def surface_angle(u):
            # work in cam's frame space
            xtangent = self.tangent_point(u)
            #print("   xtangent=", xtangent)
            rope_angle = arctan2(xtangent[1]-xpivot[1], xtangent[0]-xpivot[0])
            tangent_angle = self.tangent_angle(u)
            delta = rope_angle - tangent_angle
            #print("   rope_angle=", rope_angle, " tangent_angle=", tangent_angle);
            if delta > pi: delta = delta - 2*pi
            if delta < -pi: delta = delta + 2*pi
            return delta.A1[0]

        #print("pivot=",pivot, "  u0=", u0)
        sa0 = surface_angle(u0)
        count = 0
        step_size = 0.01
        if (sign(sa0) == 0):
            uopt=u0
        else:
            while(count < 1000):
                u1 = u0 + sign(sa0)*step_size
                sa1 = surface_angle(u1)
                #print("  u1=", u1, "  sa1=", sa1)
                if sign(sa0) != sign(sa1):
                    break
                if abs(sa0 - sa1) > pi/8:
                    step_size = step_size/4
                else:
                    u0,sa0 = u1,sa1
                    count = count + 1

            if sign(sa0) != sign(sa1):
                uopt = optimize.brentq(surface_angle, u0, u1, disp=disp)
            else:
                uopt = sa0
        if (disp):
            print ("count=", count, "uopt=", uopt, "surface_angle=", surface_angle(uopt))
        return uopt

    def tangent_angle(self, u):
        "return cam's frame-space tangent angle given cam parameter u"
        dv = splev(mod(u,1.0), self.camtck, der=1)
        angle = arctan2(dv[1], dv[0])
        return angle
    
    def integral(self, u1, u0=None):
        "return cam's definite line integral given two cam parameters"
        if u0 is None:
            u0 = 0.0*u1
        u1m = mod(u1,1.0)
        u0m = mod(u0,1.0)
        y = self.lspline(u1m) - self.lspline(u0m) + \
            (self.lspline(1.0)-self.lspline(0.0)) * ((u1-u1m) - (u0-u0m))
        if isscalar(u1):
            y = y[0]
        return y

    def eval(self, disp=False):
        #pdb.set_trace()
        self.camu = self.tangent_u(self.frame1.frame2world(self.xframe1), self.camu, disp)
        #print("camu=", self.camu)
        xtangent = self.cam_frame.frame2world(self.tangent_point(self.camu))

        self.xframe0 = self.cam_frame.world2frame(xtangent)
        self.wound_length = self.integral(1.0, self.camu)
        if (self.total_length is None):
            self.total_length = length(self.frame0.frame2world(self.xframe0) - \
                                       self.frame1.frame2world(self.xframe1)).A1[0] + \
                                self.wound_length
        self.length = self.total_length - self.wound_length  # rod length
        return (Rod.eval(self))

    def draw(self, cr, drawForces):
        self.eval(disp=True)
        cr.set_source_rgb(0.3,0.5,1.0)
        pix,foo = cr.device_to_user_distance(1.0,1.0)
        cr.set_line_width(1.0*pix)
        for u in linspace(0, 1):
            xf = self.tangent_point(u)
            xw = self.cam_frame.frame2world(xf).A1
            cr.line_to(xw[0], xw[1])
        cr.close_path()
        cr.stroke()
        xcenter  = self.cam_frame.frame2world(self.center)
        xtangent = self.cam_frame.frame2world(self.tangent_point(self.camu))
        cr.set_source_rgb(0.2,0.2,0.4)
        cr.move_to(xcenter[0], xcenter[1])
        cr.line_to(xtangent[0], xtangent[1])
        cr.stroke()
        Rod.draw(self, cr, drawForces)

def perp_distsq(a, b, c):
    "return square of distance from a to line bc"
    ba = a-b
    bc = c-b
    vec = bc * (ba.T*bc / (bc.T*bc)) - ba
    return (vec.T * vec).A1[0]
