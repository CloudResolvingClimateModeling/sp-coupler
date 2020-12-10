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
                for i,p in enumerate(points):
                    if g.contains(shapely.geometry.Point(p)): result.append(i)

                    # try also with the grid point mapped to the -180 ... 180 interval
                    # 
                    q = ((p[0] -180) % 360 - 180, p[1])
                    if g.contains(shapely.geometry.Point(q)): result.append(i)
                    
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
