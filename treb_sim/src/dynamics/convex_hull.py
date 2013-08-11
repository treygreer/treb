from numpy import arctan2, sqrt, lexsort, copy, argmin, resize, inf, vstack, shape
from pylab import find
from math import pi
import pdb

def hull(points):
    "given an (N,2) array of points, returns the convex hull array"
    hull = copy(points)
    if shape(hull)[0] <= 3:
        return hull

    # swap (first) bottom-most point to position zero
    idx = argmin(hull[:,1]) 
    hull[0],hull[idx] = hull[idx],copy(hull[0])

    idx = 0
    while(True):
        vecs = hull - hull[idx]
        angle = arctan2(vecs[:,1], vecs[:,0]) + pi
        angle[idx] = inf
        if idx>0:
            angle = angle-angle[idx-1]
            angle[idx-1] = inf
            angle[angle<0] = angle[angle<0] + 2*pi
        dist = sqrt(vecs[:,1]**2 + vecs[:,0]**2)
        nxt = find(angle == min(angle))  # major key: minimum angle
        nxt = nxt[argmin(dist[nxt])]     # minor key: minimum distance
        idx = idx+1
        if nxt == 0: break
        hull[idx],hull[nxt] = hull[nxt],copy(hull[idx])
    return (resize(hull, [idx,2]))


