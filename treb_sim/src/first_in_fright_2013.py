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

from flight import Flight

import scipy
import scipy.interpolate
import numpy as np
from math import pi, sin, cos, sqrt, acos, atan2
from scipy.optimize.optimize import fmin
from scipy.optimize.minpack import fsolve
#from scipy.interpolate.fitpack2 import UnivariateSpline
from pylab import plot

scipy.set_printoptions(precision=5, linewidth=200)
def treb( sling_length = 9.3,    # sling length, feet
          ramp_length = 11,          # ramp length, feet
          link_sum = 5.587,          # sum of upper and lower link lengths, feet
          hanger_x = 11.38508,       # feet
          hanger_y = -2,
          hinge_x = (6.+2.)/12.,     # feet
          hinge_y = -4.0,
          alpha=90,             # arm start angle, ccw from horizontal (degrees)
          omega=10,             # cocked angle between upper link and lower link
          cw_drop = 5.0,        # feet
          cw_weight = 4581 + 2000,     # pounds
          cw_moment_arm = 10.41, # distance from hinge to cw center of gravity, feet
          cw_moment = 3.516e6,   # counterweight moment about its CG, lb*ft^2
          upper_link_weight = 2*58.,  # pounds
          lower_link_weight = 2*52.,  # pounds
          link_axle_weight = 106,     # pounds
          connector_rod_weight = 84.8,  # pounds
          connector_brace_weight = 105,  # pounds
          pumpkin_weight = 10.0,     # pounds
          sling_weight = 1.7,        # pounds
          sim_duration = 2.0,  # seconds
          dry_fire = False,    # True to disable sling from time 0 
          time_step = 0.001,   # seconds
          slide_y = -9,        # feet
          arm_depth = (10.+1./4.)/12.,  # inches
          arm_thick = (5.+1./4.)/12.,   # inches
          arm_end_depth = (6.+5./8)/12.,# inches
          arm_end_thick = (3.+1./8)/12.,# inches
          release_pin_weight = 9, # pounds
          release_time = 0.0, #seconds
          debug = True):

    sim = dynamics.simulation.Simulation(max_time=sim_duration,
                                         time_step=time_step)
    sim.debug=debug

    # convert arguments to metric and radians
    sling_length = foot2meter(sling_length)
    hanger_pos = foot2meter(np.array((hanger_x, hanger_y)))
    del hanger_x, hanger_y
    hinge_pos = foot2meter(np.array((hinge_x, hinge_y)))
    del hinge_x, hinge_y
    slide_y = foot2meter(slide_y)
    arm_depth = foot2meter(arm_depth)
    arm_thick = foot2meter(arm_thick)
    arm_end_depth = foot2meter(arm_end_depth)
    arm_end_thick = foot2meter(arm_end_thick)
    ramp_length = foot2meter(ramp_length)
    link_sum = foot2meter(link_sum)
    sim.release_time = release_time
    alpha = scipy.deg2rad(alpha)
    omega = scipy.deg2rad(omega)
    cw_drop = foot2meter(cw_drop)
    cw_mass = lb2kgram(cw_weight)
    cw_moment_arm = foot2meter(cw_moment_arm)
    cw_moment = cw_moment / 32.174049 * 0.00029263965 # convert lb to slug, then
               # slug*in^2  to  kgram*meter^2
    connector_rod_mass = lb2kgram(connector_rod_weight)
    connector_brace_mass = lb2kgram(connector_brace_weight)
    upper_link_mass = lb2kgram(upper_link_weight)
    lower_link_mass = lb2kgram(lower_link_weight)
    link_axle_mass = lb2kgram(link_axle_weight)
    pumpkin_mass = lb2kgram(pumpkin_weight)
    sling_mass = lb2kgram(sling_weight)
    release_pin_mass = lb2kgram(release_pin_weight)

    # long arm length to reach slide
    long_arm_length = -slide_y / np.sin(alpha) - inch2meter(0)

    # compute rest cw position thru triangulation
    rest_cw_ctr = circle_intersection(hanger_pos, link_sum,
                                      hinge_pos, ramp_length)

    # compute cocked cw position on circle about hinge, up 'drop' meters from rest position
    cocked_cw_ctr = np.array((None, rest_cw_ctr[1] + cw_drop))
    # ramp_length**2 = (x-hinge_x)**2 + (y-hinge_y)**2
    cocked_cw_ctr[0] = hinge_pos[0] + sqrt(ramp_length**2 - (cocked_cw_ctr[1]-hinge_pos[1])**2)

    # cocked connection point is on ellipse w/ foci at hanger and cocked_cw, 'string' length
    #  equal to link_sum, 'string' interior angle omega.  In maxima:
    #    r2: s-r1
    #    eq1: d^2 = r1^2+r2^2-2*r1*r2*cos(omega)
    #    solve(eq1, r1)
    d = length_(hanger_pos - cocked_cw_ctr)
    s = link_sum
    sol1 = -(sqrt(s**2*cos(omega)**2 + 2*d**2*cos(omega)-s**2+2*d**2) - s*cos(omega) - s)/(2*cos(omega)+2)
    sol2 = (sqrt(s**2*cos(omega)**2 + 2*d**2*cos(omega)-s**2+2*d**2) + s*cos(omega) + s)/(2*cos(omega)+2)
    upper_link_length = min(sol1,sol2)
    lower_link_length = max(sol1,sol2)
    if abs((upper_link_length+lower_link_length-link_sum)/link_sum) > 0.001:
        print("link sum error")
        print("  upper_link_length=", meter2foot(upper_link_length))
        print("  lower_link_length=", meter2foot(lower_link_length))
        print("  link_sum=", meter2foot(link_sum))
        raise ValueError
    cocked_connection_pos = circle_intersection(cocked_cw_ctr, lower_link_length,
                                                hanger_pos, upper_link_length)

    # all link angles measured at top of link
    cocked_upper_link_angle = rot2radians(cocked_connection_pos - hanger_pos)
    cocked_lower_link_angle = rot2radians(cocked_cw_ctr - cocked_connection_pos)
    rest_upper_link_angle = rot2radians(rest_cw_ctr - hanger_pos)
    rest_lower_link_angle = rest_upper_link_angle
    rest_connection_pos = hanger_pos + upper_link_length * radians2rot(rest_upper_link_angle)

    # end of short arm is on ellipse with foci at axle and cocked connection, with 'string' length
    #  distance from axle to rest connection point.
    axle_rest_connection_distance = length_(rest_connection_pos)
    ellipse_axis_angle = rot2radians(-cocked_connection_pos)
    ellipse_a = axle_rest_connection_distance / 2.0
    ellipse_f = length_(cocked_connection_pos) / 2.0
    ellipse_e = ellipse_f / ellipse_a
    theta = ellipse_axis_angle - cocked_upper_link_angle
    connector_length = ellipse_a * (1-ellipse_e**2) / (1 - ellipse_e*cos(theta))

    # cocked_connection angle measured at connection point
    cocked_connection_angle = cocked_upper_link_angle
    cocked_short_arm_end = cocked_connection_pos + connector_length * radians2rot(cocked_connection_angle)
    short_arm_length = length_(cocked_short_arm_end)
    if abs((short_arm_length + connector_length - axle_rest_connection_distance)/axle_rest_connection_distance) > 0.001:
        print ("short arm length error:")
        print ("  ellipse_a=", meter2foot(ellipse_a))
        print ("  ellipse_f=", meter2foot(ellipse_f))
        print ("  ellipse_e=", ellipse_e)
        print ("  theta=", scipy.rad2deg(theta))
        print ("  connector_length=", meter2foot(connector_length))
        print ("  short_arm_length=", meter2foot(short_arm_length))
        print ("  axle_rest_connection_distance=",
                meter2foot(axle_rest_connection_distance))
        raise ValueError

    # short arm angle measured at axle
    cocked_short_arm_angle = rot2radians(cocked_short_arm_end)

    # compute beta, angle from long arm to short arm
    beta = pi + alpha - cocked_short_arm_angle

    # long arm end, cocked
    cocked_long_arm_end = long_arm_length * radians2rot(pi+alpha)

    # other dimensions
    pumpkin_diameter = inch2meter(8.0)
    pumpkin_ctr = cocked_long_arm_end + np.array((sling_length, 0.0))

    if debug:
        # rest short arm angle and position (for printing only)
        rest_short_arm_angle = rot2radians(rest_connection_pos)
        rest_short_arm_end = short_arm_length * radians2rot(rest_short_arm_angle)

        # rest long arm angle and position (for printing only)
        rest_long_arm_angle = (pi+alpha) + (rest_short_arm_angle - cocked_short_arm_angle)
        rest_long_arm_end = long_arm_length * radians2rot(rest_long_arm_angle)

        print("slide_y=", meter2foot(slide_y))
        print("long_arm_length=", meter2foot(long_arm_length))
        print("pumpkin=", meter2foot(pumpkin_ctr))
        print("hanger=", meter2foot(hanger_pos))
        print("cocked_connection=", meter2foot(cocked_connection_pos))
        print("cocked_cw=", meter2foot(cocked_cw_ctr))
        print("cocked_short_arm=", meter2foot(cocked_short_arm_end))
        print("cocked_long_arm=", meter2foot(cocked_long_arm_end))
        print("cocked_lower_link_angle=", scipy.rad2deg(cocked_lower_link_angle))
        print("rest_lower_link_angle=", scipy.rad2deg(rest_lower_link_angle))
        print("connector_length=", meter2foot(connector_length))
        print("lower_link_length=", meter2foot(lower_link_length))
        print("rest_cw_ctr=", meter2foot(rest_cw_ctr))
        print("rest_connection=", meter2foot(rest_connection_pos))
        print("rest_short_arm=", meter2foot(rest_short_arm_end))
        print("rest_long_arm=", meter2foot(rest_long_arm_end))
        
    ### Machine frame origin is at axle
    sim.machineFrame=Frame(sim, "machine", theta=0, origin=(0,0))
    sim.machineFrame.machine=Rectangle(sim.machineFrame,
                                       l=hanger_pos[0]+2.0,
                                       w=-slide_y+1.0,
                                       theta=0,
                                       origin=(hanger_pos[0]/2,
                                               (slide_y)/2),
                                       mass=lb2kgram(5000),
                                       color=(0,0,0))
    front_foot_pos = (hanger_pos[0], slide_y-0.5)
    rear_foot_pos = (0, slide_y - 0.5)
    sim.machineFrame.rear_foot=Rectangle(sim.machineFrame,
                                         l=0.3,
                                         w=0.1,
                                         origin=rear_foot_pos,
                                         mass=0,
                                         color=(0,0,0))
    sim.machineFrame.front_foot=Rectangle(sim.machineFrame,
                                         l=0.3,
                                         w=0.1,
                                         origin=front_foot_pos,
                                         mass=0,
                                         color=(0,0,0))

    ### Arm frame origin is at axle.  Framespace has long arm horizontal to the left
    sim.armFrame=Frame(sim, "arm", theta=alpha, origin=(0,0))
    sim.armFrame.long_arm=Beam(sim.armFrame,
                       x0=-long_arm_length, d0=arm_end_depth, t0=arm_end_thick,
                       x1=0,         d1=arm_depth,     t1=arm_thick,
                       density=pine_density,
                       color=(0.8,0.3,0))
    sim.armFrame.short_arm=dynamics.object.Rectangle(sim.armFrame,
                                     l=inch2meter(18.99),
                                     w=inch2meter(8.0),
                                     theta=-beta,
                                     origin=(-inch2meter(15.0)*cos(beta),
                                              inch2meter(15.0)*sin(beta)),
                                     mass=lb2kgram(53),
                                     color=(0.8,0.3,0))
    sim.armFrame.connector_pin=dynamics.object.Circle(sim.armFrame,
                                               radius=inch2meter(2.0),
                                               origin=(-short_arm_length*cos(beta),
                                                        short_arm_length*sin(beta)),
                                               mass=lb2kgram(1),
                                               color=(0.8,0.3,0))
    sim.armFrame.long_arm_plate=dynamics.object.Rectangle(sim.armFrame,
                                                   l=inch2meter(27.5),
                                                   w=inch2meter(8.0),
                                                   theta=0.0,
                                                   origin=(inch2meter(-6.25), 0),
                                                   mass=lb2kgram(63),
                                     color=(0.8,0.3,0))
    sim.armFrame.release_pin=dynamics.object.Circle(sim.armFrame,
                                             radius=inch2meter(6),
                                             origin=(-long_arm_length, 0),
                                             mass=release_pin_mass, color=(1.0, 1.0, 1.0))

    # Wdight frame origin is at pivot point, ramp horizontal to the right
    cocked_ramp_angle = rot2radians(cocked_cw_ctr-hinge_pos)
    sim.weightFrame=dynamics.frame.Frame(sim, "weight", theta=cocked_ramp_angle, origin=hinge_pos)
    sim.weightFrame.ramp = dynamics.object.Rectangle(sim.weightFrame, l=ramp_length, w=inch2meter(4),
                                            mass=0, color=(0.3,0.5,0.2),
                                            origin = (ramp_length/2,0))
    sim.weightFrame.cw = dynamics.object.Rectangle(sim.weightFrame, l=foot2meter(2.6), w=foot2meter(2.6),
                                            color=(0.3,0.5,0.2),
                                            mass=cw_mass,
                                            origin = (cw_moment_arm,0),
                                            moment = cw_moment)

    # Lower link frame origin is at end of ramp
    sim.lowerLinkFrame = dynamics.frame.Frame(sim, "lower link", origin=cocked_cw_ctr,
                                        theta = cocked_lower_link_angle-pi)
    sim.lowerLinkFrame.link = dynamics.object.Rectangle(sim.lowerLinkFrame, l=lower_link_length, w=inch2meter(6),
                                                 mass=lower_link_mass, color=(1.0,0.0,0.0),
                                                 origin=(lower_link_length/2, 0.0))
    sim.lowerLinkFrame.axle=dynamics.object.Circle(sim.lowerLinkFrame,
                                            radius=inch2meter(3),
                                            origin=(lower_link_length, 0.0),
                                            mass=link_axle_mass, color=(1.0, 0.0, 0.0))

    # Upper link frame origin is the hanger
    cocked_upper_link_angle = rot2radians(cocked_connection_pos-hanger_pos)
    sim.upperLinkFrame = dynamics.frame.Frame(sim, "upper link", origin=hanger_pos,
                                        theta = cocked_upper_link_angle)
    sim.upperLinkFrame.link = dynamics.object.Rectangle(sim.upperLinkFrame, l=upper_link_length, w=inch2meter(6),
                                                 mass=upper_link_mass, color=(1.0,0.0,0.0),
                                                 origin=(upper_link_length/2, 0.0))

    # Connector frame origin is the end of the short arm
    sim.connectorFrame = dynamics.frame.Frame(sim, "connector", origin=cocked_short_arm_end,
                                        theta = rot2radians(cocked_connection_pos - cocked_short_arm_end))
    sim.connectorFrame.rod = dynamics.object.Rectangle(sim.connectorFrame, l=connector_length,
                                                      w=inch2meter(2),
                                                      mass=connector_rod_mass,
                                                      color=(0.0, 0.0, 0.0),
                                                      origin=(connector_length/2, 0.0))
    sim.connectorFrame.stiffener = dynamics.object.Rectangle(sim.connectorFrame, l=connector_length,
                                                      w=inch2meter(4.0),
                                                      mass=lb2kgram(100),
                                                      color=(0.0, 0.0, 0.0),
                                                      origin=(connector_length/2, inch2meter(3.0)))
    sim.connectorFrame.brace = dynamics.object.Rectangle(sim.connectorFrame, l=foot2meter(2),
                                                       w=inch2meter(4),
                                                       mass=connector_brace_mass,
                                                       color=(0.0, 0.0, 0.0),
                                                       origin=(connector_length-foot2meter(1), 0.0))

    # Pumpkin
    sim.pumpkinFrame=dynamics.frame.Frame(sim, "pumpkin", origin=pumpkin_ctr)
    sim.pumpkinFrame.pumpkin=dynamics.object.Circle(sim.pumpkinFrame,
                                             radius=pumpkin_diameter/2.0,
                                             mass=pumpkin_mass, color=(1.0, 0.5, 0))
    sim.pumpkinFrame.sling=dynamics.object.Circle(sim.pumpkinFrame,
                                             radius=pumpkin_diameter/2.0,
                                             mass=sling_mass, color=(1.0, 0.5, 0))

    # initialize frames
    for frame in sim.frames:
        frame.init()

    # define constraints
    sim.rear_foot = Nail(sim, "rear foot",
                         obj=sim.machineFrame.rear_foot,
                         xobj=(0,0),
                         xworld=rear_foot_pos)

    sim.front_foot = NailSpring(sim, "front foot",
                         obj=sim.machineFrame.front_foot,
                         xobj=(0,0),
                         x_world=front_foot_pos,
                         spring_constant=1e6,
                         damping_constant=500e3)

    sim.axle = Pin(sim, "axle",
                             obj0=sim.armFrame.long_arm,
                             xobj0=(0, 0),
                             obj1=sim.machineFrame)

    sim.hinge =Pin(sim, "hinge",
                              obj0=sim.weightFrame.ramp,
                              xobj0=(-ramp_length/2, 0.0),
                              obj1=sim.machineFrame)

    sim.hanger = Pin(sim, "hanger",
                               obj0=sim.upperLinkFrame.link,
                               xobj0=(-upper_link_length/2.0,0.0),
                               obj1=sim.machineFrame)

    sim.linkPin = Pin(sim, "linkPin",
                               obj0=sim.upperLinkFrame.link,
                               xobj0= (upper_link_length/2.0, 0.0),
                               obj1=sim.lowerLinkFrame.link,
                               xobj1 = (lower_link_length/2.0, 0.0))

    sim.rampPin = dynamics.constraint.Pin(sim, "rampPin",
                               obj0=sim.weightFrame.ramp,
                               xobj0= (ramp_length/2.0, 0.0),
                               obj1=sim.lowerLinkFrame.link,
                               xobj1 = (-lower_link_length/2.0, 0.0))

    sim.connectorPin1 = Pin(sim, "connectorPin1",
                                     obj0=sim.armFrame.connector_pin,
                                     xobj0=(0.0,0.0),
                                     obj1=sim.connectorFrame.rod,
                                     xobj1 = (-connector_length/2.0, 0.0))

    sim.connectorPin2 = Pin(sim, "connectorPin2",
                                     obj0=sim.upperLinkFrame.link,
                                     xobj0=(upper_link_length/2.0,0.0),
                                     obj1=sim.connectorFrame.rod,
                                     xobj1 = (connector_length/2.0, 0.0))

    sim.sling=Rod(sim, "sling",
                           obj0=sim.armFrame.long_arm, xobj0=(-long_arm_length,
                                                              0),
                           obj1=sim.pumpkinFrame.pumpkin, xobj1=(0.0,0.0),
                           length=sling_length)
    '''
    sim.trigger = Rod(sim, "trigger",
                               obj0=sim.pumpkinFrame.pumpkin,
                               xobj0= (0.0, 0.0),
                               obj1=sim.machineFrame.front_foot,
                               xobj1= (0.0,0.0))
    '''

    sim.slide=Shelf(sim, "slide",
                             obj=sim.pumpkinFrame.pumpkin,
                             xobj=(0,0),
                             height=slide_y)

    if (dry_fire):
        sim.sling.enabled = False

    print( "    running simulation")
    from time import clock
    tstart=clock()
    sim.run(continue_sim, debug=debug)
    print ("    done: time=%g sec" % (clock()-tstart))

    if not sim.release_time:
        sim.range = Y2range(sim,sim.Y)
        range_spline = scipy.interpolate.UnivariateSpline(sim.t, sim.range, k=3,s=0.0)
        d0,t0 = max( (range,time) for range,time in zip(sim.range, sim.t) ) # find guess
        sim.tmax = fsolve(range_spline, t0, args=1)  # root of first derivative of range
        sim.maxrange = range_spline(sim.tmax)
        launchDegrees_spline = scipy.interpolate.UnivariateSpline(sim.t, Y2launchDegrees(sim.Y), k=3,s=0.0)
        sim.launchDegrees = launchDegrees_spline(sim.tmax)
        print ("     distance=%g feet at %g sec" % (meter2foot(sim.maxrange), sim.tmax))
    else:
        sim.range=np.zeros(len(sim.t))
        sim.maxrange=0

    sim.Fmax = max(sim.hanger.Fvec())
    print("     max force on hanger = %g pounds" % (newton2lb(sim.Fmax)))
    return(sim)

def circle_intersection(ctr1, rad1, ctr2, rad2):
    """Return intersection of two circles.

    Intersection returned is the one in the ccw direction from the vector
    ctr1->ctr2.

    """

    base_len = length_(ctr2-ctr1)
    # alpha is angle from vector ctr1->ctr2 to vector ctr1->isect
    alpha = acos( (base_len**2 + rad1**2 - rad2**2) / (2 * base_len * rad1) )
    # beta is angle from positive x axis to vector ctr1->ctr2
    beta = rot2radians(ctr2-ctr1)
    isect = ctr1 + rad1*radians2rot(alpha+beta)
    return isect

def continue_sim(sim, time, y):
    "continue simulation?"

    #if time>0.001:
    #    sim.trigger.enabled = False

    if sim.slide.enabled:
        shelf_force = sim.slide.forces[0][1]
        if shelf_force < 0.0:
            sim.slide.enabled = False

    if 0:
        if sim.sling.enabled:
            v = sim.pumpkinFrame.v
            angle = atan2(v.A[1], v.A[0])
            if v.A[0] > 0.0 and v.A[1] > 0.0 and angle <= sim.release_angle:
                sim.maxrange = Y2range(sim,y)[0]
                sim.sling.enabled = False
                #return False
        return True
    else:
        if sim.release_time:
            if time >= sim.release_time:
                sim.sling.enabled = False
            return True
        if sim.armFrame.theta >= -3*pi/4:
            return True
        if sim.pumpkinFrame.v.A1[1] > 0:
            return True
        return False

def Y2range(sim, Y, with_air_friction=True):
    if (len(Y.shape)==1):
        Y = Y.reshape([1,len(Y)])
    idx = sim.pumpkinFrame.idx
    x0 = Y[:,6*idx]
    y0 = Y[:,6*idx+1]
    vx0 = Y[:,6*idx+3]
    vy0 = Y[:,6*idx+4]
    
    if not with_air_friction:
        tof = 2.0 * vy0 / scipy.constants.g
        tof[tof<0.0] = 0.0
        return (tof*vx0)
    else:        
        range = np.zeros(len(x0))
        flight = Flight(mass=sim.pumpkinFrame.pumpkin.mass,
                        area=pi*sim.pumpkinFrame.pumpkin.radius**2)
        for i in np.arange(len(x0)):
            if (vy0[i] > 0) & (vx0[i] > 0):
                flight.run([x0[i],y0[i]], [vx0[i],vy0[i]])
                range[i] = flight.range()
        return range

def Y2launchDegrees(Y):
    if (len(Y.shape)==1):
        Y = Y.reshape([1,len(Y)])
    vx = Y[:,33]
    vy = Y[:,34]
    return (180./pi*np.arctan2(vy, vx))


def trebPEvec(sim):
    return (sim.machineFrame.PEvec()+
            sim.weightFrame.PEvec() +
            sim.upperLinkFrame.PEvec() +
            sim.lowerLinkFrame.PEvec() +
            sim.connectorFrame.PEvec() +
            sim.armFrame.PEvec())

def trebKEvec(sim):
    return (sim.machineFrame.KEvec() +
            sim.weightFrame.KEvec() +
            sim.upperLinkFrame.KEvec() +
            sim.lowerLinkFrame.KEvec() +
            sim.connectorFrame.KEvec() +
            sim.armFrame.KEvec())

def plotEnergies(sim):
    plot (sim.t, trebPEvec(sim) - min(trebPEvec(sim)))
    plot (sim.t, trebKEvec(sim))
    plot (sim.t, (trebPEvec(sim) - min(trebPEvec(sim)) +
                  trebKEvec(sim)))
    plot (sim.t, trebKEvec(sim))
    plot (sim.t, sim.pumpkinFrame.KEvec() + sim.pumpkinFrame.PEvec())
    

def opt(X):
    global sim, X0
    X0 = X 
    print ("X=", X)
    try:
        sim = treb(debug=False, time_step=0.0001, sim_duration=0.8,
                   sling_length=X[0])
        return -sim.maxrange
        #return -sim.maxrange / sim.Fmax**0.10

    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except:
        return 0.0

#X0 =  array([  8.70381,   6.08564,  10.3123 ])
#X0 =  array([  8,   6,  10 ])
#X0 = [ 9.62859,  6.23794,  9.98966]
#X0 = [  8.70153,   6.04452,  10.43426]
#X0 = array([  8.68625,   6.00475,  10.44   ])
#X0 = array([  8.21222,   5.58682,  11.43518, -9.0])
#X0 = array([8.411, 5.587, 11.433])
X0 = np.array([9.3])
#lower =  array([ 6.0,      3.0,      5.0])
#upper =  array([ 12.0,     9.0,      12.0])
#result=scipy.optimize.fmin(opt, X0)
#result=scipy.optimize.fmin_l_bfgs_b(opt, X0, approx_grad=True, bounds=None)
#result=scipy.optimize.anneal(opt, X0, lower=lower, upper=upper, T0=0.001, feps=1e-60, full_output=True)

if __name__ == '__main__':
    sim=treb(debug=True)
    anim=Animation(sim, Y2range)