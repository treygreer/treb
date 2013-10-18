''' a bit more in the comment...
'''

import dynamics.simulation
from dynamics.frame import Frame
from dynamics.object import Rectangle, Circle, Beam
from dynamics.constraint import Nail, Rod, Pin, Shelf

from dynamics.constants import foot2meter, inch2meter, meter2foot
from dynamics.misc import length_, rot2radians, radians2rot
from dynamics.constants import lb2kgram, kgram2lb, newton2lb
from dynamics.constants import pine_density, steel_density

import scipy
from scipy import constants
from numpy import pi, sqrt, array, sin, cos, arccos, arctan2
from scipy import deg2rad, rad2deg
#from scipy.optimize.optimize import fmin
#from scipy.optimize.minpack import fsolve
#from scipy.interpolate.fitpack2 import UnivariateSpline
from pylab import plot

scipy.set_printoptions(precision=5, linewidth=200)
def treb( sling_length = 8.54665,              # sling length, feet
          ramp_length = 11,                    # ramp length, feet
          link_sum = 5.587,                    # sum of upper and lower link lengths, feet
          hanger_x = 11.38508,       # feet
          hanger_y = -2,
          hinge_x = (6.+2.)/12.,    # feet
          hinge_y = -4.0,
          alpha=90,              # arm start angle, ccw from horizontal (degrees)
          omega=10,              # cocked angle between upper link and lower link (degrees)
          cw_drop = 5.0,         # feet
          cw_weight = 4000.0,    # pounds
          long_arm_weight9 = 70,  # weight at nominal 9 foot length, pounds
          short_arm_weight = 70.0,   # pounds
          upper_link_in2 = 2*6.0,  # cross section, sq inches
          lower_link_in2 = 2*6.0,  # cross section, sq inches
          connector_in2 = pi*(2.0/2.0)**2, # cross section, sq inches
          ramp_weight = 400, # pounds
          pumpkin_weight = 10.0,     # pounds
          sim_duration = 2.0,  # seconds
          release_angle = 40,  # pumpkin release angle, degrees above horizon
          dry_fire = False,    # True to disable sling from time 0
          time_step = 0.001,   # seconds
          slide_y = -9,        # feet
          arm_depth = (10.+1./4.)/12.,
          arm_thick = (5.+1./4.)/12.,
          arm_end_depth = (6.+5./8)/12.,
          arm_end_thick = (3.+1./8)/12.,
          debug = True):

    sim = dynamics.simulation.Simulation(max_time=sim_duration, time_step=time_step)
    sim.debug=debug

    # convert arguments to metric and radians
    sling_length = foot2meter(sling_length)
    hanger_pos = foot2meter(array((hanger_x, hanger_y)))
    del hanger_x, hanger_y
    hinge_pos = foot2meter(array((hinge_x, hinge_y)))
    del hinge_x, hinge_y
    slide_y = foot2meter(slide_y)
    arm_depth = foot2meter(arm_depth)
    arm_thick = foot2meter(arm_thick)
    arm_end_depth = foot2meter(arm_end_depth)
    arm_end_thick = foot2meter(arm_end_thick)
    ramp_length = foot2meter(ramp_length)
    link_sum = foot2meter(link_sum)
    sim.release_angle = deg2rad(release_angle)
    del release_angle
    alpha = deg2rad(alpha)
    omega = deg2rad(omega)
    cw_drop = foot2meter(cw_drop)
    cw_mass = lb2kgram(cw_weight)
    ramp_mass = lb2kgram(ramp_weight)
    #long_arm_mass9 = lb2kgram(long_arm_weight9)
    short_arm_mass = lb2kgram(short_arm_weight)
    connector_m2 = connector_in2 * inch2meter(1)**2
    upper_link_m2 = upper_link_in2 * inch2meter(1)**2
    lower_link_m2 = lower_link_in2 * inch2meter(1)**2
    pumpkin_mass = lb2kgram(pumpkin_weight)

    # long arm length to reach slide
    long_arm_length = -slide_y / sin(alpha)
    #long_arm_mass = long_arm_mass9 * meter2foot(long_arm_length)/9

    # compute rest cw position thru triangulation
    rest_cw_ctr = circle_intersection(hanger_pos, link_sum,
                                      hinge_pos, ramp_length)

    # compute cocked cw position on circle about hinge, up 'drop' meters from rest position
    cocked_cw_ctr = array((None, rest_cw_ctr[1] + cw_drop))
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
        print ("  theta=", rad2deg(theta))
        print ("  connector_length=", meter2foot(connector_length))
        print ("  short_arm_length=", meter2foot(short_arm_length))
        print ("  axle_rest_connection_distance=", meter2foot(axle_rest_connection_distance))
        raise ValueError

    # short arm angle measured at axle
    cocked_short_arm_angle = rot2radians(cocked_short_arm_end)

    # compute beta, angle from long arm to short arm
    beta = pi + alpha - cocked_short_arm_angle

    # long arm end, cocked
    cocked_long_arm_end = long_arm_length * radians2rot(pi+alpha)

    # other dimensions
    pumpkin_diameter = inch2meter(8.0)
    pumpkin_ctr = cocked_long_arm_end + array((sling_length, 0.0))

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
        print("cocked_lower_link_angle=", rad2deg(cocked_lower_link_angle))
        print("rest_lower_link_angle=", rad2deg(rest_lower_link_angle))
        print("connector_length=", meter2foot(connector_length))
        print("lower_link_length=", meter2foot(lower_link_length))
        print("rest_cw_ctr=", meter2foot(rest_cw_ctr))
        print("rest_connection=", meter2foot(rest_connection_pos))
        print("rest_short_arm=", meter2foot(rest_short_arm_end))
        print("rest_long_arm=", meter2foot(rest_long_arm_end))

    ### Arm frame origin is at axle.  Framespace has long arm horizontal to the left
    sim.armFrame=Frame(sim, "arm", theta=alpha, origin=(0,0))
    sim.armFrame.long_arm=Beam(sim.armFrame,
                       x0=-long_arm_length, d0=arm_end_depth, t0=arm_end_thick,
                       x1=0,         d1=arm_depth,     t1=arm_thick,
                       density=pine_density,
                       color=(0.8,0.3,0))
    sim.armFrame.short_arm=Rectangle(sim.armFrame,
                                     l=short_arm_length,
                                     w=inch2meter(8),
                                     theta=-beta,
                                     origin=(-short_arm_length/2*cos(beta),
                                             short_arm_length/2*sin(beta)),
                                     mass=short_arm_mass,
                                     color=(0.8,0.3,0))

    # Ramp frame origin is at pivot point, ramp horizontal to the right
    cocked_ramp_angle = rot2radians(cocked_cw_ctr-hinge_pos)
    sim.rampFrame=Frame(sim, "ramp", theta=cocked_ramp_angle, origin=hinge_pos)
    sim.rampFrame.ramp = Rectangle(sim.rampFrame, l=ramp_length, w=inch2meter(4),
                                            mass=ramp_mass, color=(0.3,0.5,0.2),
                                            origin = (ramp_length/2,0))
    sim.rampFrame.cw = Rectangle(sim.rampFrame, l=foot2meter(2.6), w=foot2meter(2.6),
                                   mass=cw_mass, color=(0.3,0.5,0.2),
                                   origin = (ramp_length,0))

    # Lower link frame origin is at end of ramp
    lower_link_mass = lower_link_length * lower_link_m2 * steel_density
    sim.lowerLinkFrame = Frame(sim, "lowerLink", origin=cocked_cw_ctr,
                                        theta = cocked_lower_link_angle-pi)
    sim.lowerLinkFrame.link = Rectangle(sim.lowerLinkFrame, l=lower_link_length, w=inch2meter(6),
                                                 mass=lower_link_mass, color=(1.0,0.0,0.0),
                                                 origin=(lower_link_length/2, 0.0))

    # Upper link frame origin is the hanger
    cocked_upper_link_angle = rot2radians(cocked_connection_pos-hanger_pos)
    upper_link_mass = upper_link_length * upper_link_m2 * steel_density
    sim.upperLinkFrame = Frame(sim, "upperLink", origin=hanger_pos,
                                        theta = cocked_upper_link_angle)
    sim.upperLinkFrame.link = Rectangle(sim.upperLinkFrame, l=upper_link_length, w=inch2meter(6),
                                                 mass=upper_link_mass, color=(1.0,0.0,0.0),
                                                 origin=(upper_link_length/2, 0.0))

    # Connector frame origin is the end of the short arm
    connector_mass = connector_length * connector_m2 * steel_density
    sim.connectorFrame = Frame(sim, "connector", origin=cocked_short_arm_end,
                                        theta = rot2radians(cocked_connection_pos - cocked_short_arm_end))
    sim.connectorFrame.connector = Rectangle(sim.connectorFrame, l=connector_length,
                                                      w=inch2meter(2),
                                                      mass=connector_mass,
                                                      color=(0.0, 0.0, 0.0),
                                                      origin=(connector_length/2, 0.0))
    sim.connectorFrame.crossbeam = Rectangle(sim.connectorFrame,
                                                      l=connector_length/3.0, w=inch2meter(4),
                                                      mass=lb2kgram(200), color=(0.0, 0.0, 0.0),
                                                      origin=(connector_length*5./6., 0.0))


    if (debug):
        print("connector weight = %g pounds" % kgram2lb(connector_mass))
        print("lower link weight = %g pounds" % kgram2lb(lower_link_mass))
        print("upper link weight = %g pounds" % kgram2lb(upper_link_mass))
    
    # Pumpkin
    sim.pumpkinFrame=Frame(sim, "pumpkin", origin=pumpkin_ctr)
    sim.pumpkinFrame.pumpkin=Circle(sim.pumpkinFrame,
                                             radius=pumpkin_diameter/2.0,
                                             mass=pumpkin_mass, color=(1.0, 0.5, 0))

    # initialize frames
    for frame in sim.frames:
        frame.init()

    # define constraints
    sim.axle = Nail(sim, "axle",
                             obj=sim.armFrame.long_arm,
                             xobj=(0, 0),
                             xworld=(0.0, 0.0))

    sim.hinge =Nail(sim, "hinge",
                              obj=sim.rampFrame.ramp,
                              xobj=(-ramp_length/2, 0.0),
                              xworld=hinge_pos)

    sim.hanger = Nail(sim, "hanger",
                               obj=sim.upperLinkFrame.link,
                               xobj=(-upper_link_length/2.0,0.0),
                               xworld=hanger_pos)

    sim.linkPin = Pin(sim, "linkPin",
                               obj0=sim.upperLinkFrame.link,
                               xobj0= (upper_link_length/2.0, 0.0),
                               obj1=sim.lowerLinkFrame.link,
                               xobj1 = (lower_link_length/2.0, 0.0))

    sim.rampPin = Pin(sim, "rampPin",
                               obj0=sim.rampFrame.ramp,
                               xobj0= (ramp_length/2.0, 0.0),
                               obj1=sim.lowerLinkFrame.link,
                               xobj1 = (-lower_link_length/2.0, 0.0))

    sim.connectorPin1 = Pin(sim, "connectorPin1",
                                     obj0=sim.armFrame.short_arm,
                                     xobj0=(-short_arm_length/2.0,0.0),
                                     obj1=sim.connectorFrame.connector,
                                     xobj1 = (-connector_length/2.0, 0.0))

    sim.connectorPin2 = Pin(sim, "connectorPin2",
                                     obj0=sim.upperLinkFrame.link,
                                     xobj0=(upper_link_length/2.0,0.0),
                                     obj1=sim.connectorFrame.connector,
                                     xobj1 = (connector_length/2.0, 0.0))

    sim.sling=Rod(sim, "sling",
                           obj0=sim.armFrame.long_arm, xobj0=(-long_arm_length, 0),
                           obj1=sim.pumpkinFrame.pumpkin, xobj1=(0.0,0.0),
                           length=sling_length)

    sim.slide=Shelf(sim, "slide",
                             obj=sim.pumpkinFrame.pumpkin,
                             xobj=(0,0),
                             height=slide_y)

    if (dry_fire):
        sim.sling.enabled = False

    print( "    running simulation")
    from time import clock
    tstart=clock()
    sim.distance=0
    sim.trelease=0
    sim.run(continue_sim, debug=debug)
    print ("    done: time=%g sec" % (clock()-tstart))

    #dist_spline = scipy.interpolate.UnivariateSpline(sim.t, dist(sim.Y), k=3,s=0.0)
    #d0,t0 = max( (dist,time) for time,dist in zip(sim.t, dist(sim.Y)))
    #sim.tmax = fsolve(dist_spline, t0, args=1.0)  # root of first derivative of range
    #sim.maxdist = dist_spline(sim.tmax)

    print("     distance=%g feet at %g sec" % (meter2foot(sim.distance), sim.trelease))

    sim.Fmax = max(sim.hanger.Fvec())
    print("     max force on hanger = %g pounds" % (newton2lb(sim.Fmax)))

    #sim.fall = weight_origin[1] + (short_arm_length + sim.weight_arm)
    #sim.efficiency = sim.maxdist / (2*sim.fall*sim.weightFrame.mass/sim.pumpkinFrame.mass)
    return(sim)

def circle_intersection(ctr1, rad1, ctr2, rad2):
    """Return intersection of two circles.

    Intersection returned is the one in the ccw direction from the vector 
    ctr1->ctr2.

    """

    base_len = length_(ctr2-ctr1)
    # alpha is angle from vector ctr1->ctr2 to vector ctr1->isect
    alpha = arccos( (base_len**2 + rad1**2 - rad2**2) / (2 * base_len * rad1) )
    # beta is angle from positive x axis to vector ctr1->ctr2
    beta = rot2radians(ctr2-ctr1)
    isect = ctr1 + rad1*radians2rot(alpha+beta)
    return isect
        
def continue_sim(sim, t, y):
    "continue simulation?"

    if sim.slide.enabled:
        shelf_force = sim.slide.forces[0][1]
        if shelf_force < 0.0:
            sim.slide.enabled = False

    if sim.sling.enabled:
        v = sim.pumpkinFrame.v
        angle = arctan2(v.A[1], v.A[0])
        if v.A[0] > 0.0 and v.A[1] > 0.0 and angle <= sim.release_angle:
            sim.distance = dist(y)
            sim.trelease = t
            sim.sling.enabled = False
            #return False
    return True


#    if sim.armFrame.theta >= -3*pi/4:
#        return True
#    if sim.pumpkinFrame.v.A1[1] > 0:
#        return True
#    return False

def dist(Y):
    if (len(Y.shape)==1):
        Y = Y.reshape([1,len(Y)])
    vx = Y[:,33]
    vy = Y[:,34]
    tof = 2.0 * vy / constants.g
    tof[tof<0.0] = 0.0
    return (tof*vx)


def trebPEvec(sim):
    return (sim.rampFrame.PEvec() +
            sim.upperLinkFrame.PEvec() +
            sim.lowerLinkFrame.PEvec() +
            sim.connectorFrame.PEvec() +
            sim.armFrame.PEvec())

def trebKEvec(sim):
    return (sim.rampFrame.KEvec() +
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
    print( "X=", X)
    try:
        sim = treb(debug=False, time_step=0.0001, sim_duration=0.7,
                   sling_length=X[0], link_sum=X[1], hanger_x=X[2], slide_y=-9)
        #return -sim.distance
        return -sim.distance / sim.Fmax**0.10

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
X0 = array([8.54665,   5.587,    11.38508])
#lower =  array([ 6.0,      3.0,      5.0])
#upper =  array([ 12.0,     9.0,      12.0])
#result=scipy.optimize.fmin(opt, X0)
#result=scipy.optimize.fmin_l_bfgs_b(opt, X0, approx_grad=True, bounds=None)
#result=scipy.optimize.anneal(opt, X0, lower=lower, upper=upper, T0=0.001, feps=1e-60, full_output=True)

if __name__ == '__main__':
    sim=treb(debug=True)
    anim=dynamics.animation.Animation(sim, dist)