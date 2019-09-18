import numpy
import os
from omuse.units import units
from amuse.units.quantities import to_quantity

import logging
from . import haversine
import shapely.geometry

# Logger
log = logging.getLogger(__name__)

# Physical constants
pref0 = 1e5 | units.Pa    # Pa reference pressure
rd    = 287.04 | units.J/units.kg/units.K # gas constant for dry air.  J/kg K.
rv    = 461.5 | units.J/units.kg/units.K   # gas constant for water vapor. J/kg K.
cp    = 1004. | units.J/units.kg/units.K   # specific heat at constant pressure (dry air).
rlv   = 2.53e6 | units.J/units.kg   # latent heat for vaporisation
grav  = 9.81  | units.m/units.s**2  # gravity acceleration. m/s^2
mair  = 28.967 | units.g/units.mol # Molar mass of air

# Root mean square
def rms(a):
    return numpy.sqrt(numpy.mean(a ** 2))


# Exner function
def exner(p):
    return (p / pref0) ** (rd / cp)


# inverse Exner function
def iexner(p):
    return (p / pref0) ** (-rd / cp)


# sort the indices of points by distance to the target
# points and target are given in geographic coordinates (lat, lon)
# for the n closest points use find_closest_points(points, target)[:n]
def find_closest_points(points, target):
    dists = [haversine.haversine(p, target) for p in points]
    return numpy.argsort(dists)


# Retrieves the super-parametrized indices from the input mask geometries.
def get_mask_indices(points, mask_geoms, nmax=-1):
    
    if nmax == 0:
        return []  # requested no points
    result = []
    # if a single point is specified, select the Nmax closest grid points
    if len(mask_geoms) == 1 and isinstance(mask_geoms[0], shapely.geometry.Point):
        g = mask_geoms[0]
        dists = [haversine.haversine((p[0], p[1]), (g.x, g.y)) for p in points]
        return numpy.argsort(dists)[:nmax] if nmax > 0 else [numpy.argsort(dists)[0]]
    else:
        # many geometries or not a single point. points now select only one grid index.
        for g in mask_geoms:
            if isinstance(g, shapely.geometry.Point):
                dists = [haversine.haversine((p[0], p[1]), (g.x, g.y)) for p in points]
                result.append(numpy.argmin(dists))
            else:
                for i in range(len(points)):
                    if g.contains(shapely.geometry.Point(points[i])): result.append(i)

                    # try also with the grid point mapped to the -180 ... 180 interval
                    # 
                    p = ((points[i][0] -180) % 360 - 180, points[i][1])
                    if g.contains(shapely.geometry.Point(p)): result.append(i)
                    
        result = list(set(result))  # remove duplicates
        return result


# Links contents of input directory
def link_dir(files, workdir):
    os.makedirs(workdir)
    for f in files:
        fname = os.path.basename(f)
        os.symlink(os.path.abspath(f), os.path.join(workdir, fname))

def interp(x,xp,fp,**kwargs):
    x=to_quantity(x)
    xp=to_quantity(xp)
    fp=to_quantity(fp)
    return numpy.interp(x.value_in(x.unit), xp.value_in(x.unit), fp.number, **kwargs) | fp.unit
  
def searchsorted(a,v,**kwargs):
    a=to_quantity(a)
    v=to_quantity(v)
    return numpy.searchsorted(a.value_in(a.unit), v.value_in(a.unit), **kwargs)


def integral (a, b, z, q, w=None):
    """
    Calculate the integral from a to b of the piece-wise constant function q(z).
    q(z) has the value q[i] on the interval from z[i] to z[i+1].
    
    Appropriate for integrating finite-volume quantities.
    For DALES arrays, z are half-level heights (zf).


    Parameters
    ----------
    a, b : interval end points
    z : increasing array of point coordinates
    q : array of function values (length one less than z)
    w : weight array, shaped as q. Optional.
    
    """
    if len(z) != len(q) + 1:
        print("len(z) should be len(q) + 1. len(z)=%d, len(q) = %d",(len(z), len(q)))
    if a < z[0] or a > z[-1] or b < z[0] or b > z[-1]:
        print("integral: Interval end point outside range.")
        return None

    sign = 1
    if a > b:
        sign = -1
        a,b = b,a

    ia = 0    # begin index
    while z[ia+1] < a:
        ia += 1
    ib = ia   # end index. ib >= ia
    while z[ib+1] < b:
        ib += 1

    # z[ia] <= a <= z[ia+1]   and   a <= b
    # z[ib] <= b <= z[ib+1]
    
        
    #print("integral: a=%f, ia=%d, z[ia],z[ia+1] = %f, %f"%(a, ia, z[ia], z[ia+1]))
    #print("integral: b=%f, ib=%d, z[ib],z[ib+1] = %f, %f"%(b, ib, z[ib], z[ib+1]))


    # sum intervals, including the full edge intervals
    #S = 0
    #for i in range(ia, ib+1):
    #    S += q[i] * (z[i+1]-z[i])

    if w is None:
        # numpy version - sum intervals, including the full edge intervals
        S = (q[ia:ib+1] * (z[ia+1:ib+2] - z[ia:ib+1])).sum() 
        
        # subtract edge intervals    
        Sa = q[ia] * (a - z[ia])
        Sb = q[ib] * (z[ib+1] - b)

        return (S - Sa - Sb) * sign
    else:
        S = (w[ia:ib+1] * q[ia:ib+1] * (z[ia+1:ib+2] - z[ia:ib+1])).sum()         
        # subtract edge intervals    
        Sa = w[ia] * q[ia] * (a - z[ia])
        Sb = w[ib] * q[ib] * (z[ib+1] - b)

        Sw = (w[ia:ib+1] * (z[ia+1:ib+2] - z[ia:ib+1])).sum()         
        # subtract edge intervals    
        Swa = w[ia] * (a - z[ia])
        Swb = w[ib] * (z[ib+1] - b)
        return  (S - Sa - Sb) / (Sw - Swa - Swb) * sign
        

# Optimizations:
# * construct matrix of weights (sparse or with start, end indices)
#   weights depend on rho, Zh, zh
#   same weights can be used for all quantities
#
# * use end index of previous integration as starting index for next
#
# * when searching for end index, start at begin index DONE

def interp_c(Zh, zh, q, rho):
    
    """ conservative interpolation from fine to coarse levels.
    Zh is in descending order, zh in ascending order.
    returns array Q: Q[i] = integral of q(z) from Zh[i+1] to Zh[i] weighted by rho(z)

    """
    zz = zh.value_in(units.m) # dropping units for speed
    ZZ = Zh.value_in(units.m)
    qq = q.number
    rho2 = rho.number
    
    Q = numpy.zeros(len(Zh)-1)
    for i in range(len(Q)):
        if Zh[i] < zh[-1]:        
            Q[i] = integral (ZZ[i+1], ZZ[i], zz, qq, rho2)
    return Q | q.unit

def interp_rho(Zh, zh, rho):
    """ interpolate a density rho to a coarser grid """
    RHO = numpy.zeros(len(Zh)-1) | rho.unit
    for i in range(len(RHO)):
        if Zh[i] < zh[-1]:
            RHO[i] = integral (Zh[i+1], Zh[i], zh, rho) / (Zh[i] - Zh[i+1])
    return RHO
