from numpy import mat, sqrt, array, hypot, asarray
import math 

def radians2rot(theta):
    """Converts an angle [radians] to array [cos, sin]."""
    return (array([math.cos(theta), math.sin(theta)]).T)

def rot2radians(r):
    """Returns the angle [radians] of a 2-element vector."""
    return (math.atan2(r[1], r[0]))

def rot2matrix(r):
    """Converts a rotation vector to a rotation matrix."""
    r = r / hypot(r[0],r[1])  # normalize
    return mat([[r[0], -r[1]],
                [r[1],  r[0]]])
    
def rotprod(a, b):
    "return product of two rotations a and b"
    # a rotation is like a complex a+bj, so multiply like two complex quantities:
    return(array([a[0]*b[0] - a[1]*b[1],
                  a[0]*b[1] + a[1]*b[0]]))

def dual(x):
    "return 2x2 matrix dual of scalar x"
    return(mat([[0, -x],
                [x,  0]]))

def length_(x):
    "return length of a vector"
    x = asarray(x)
    return sqrt(sum(x*x))
