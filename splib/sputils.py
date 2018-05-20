import numpy
import os
import logging
import haversine
import shapely.geometry

# Logger
log = logging.getLogger(__name__)

# Physical constants
pref0 = 1e5  # Pa reference pressure
rd = 287.04  # gas constant for dry air.  J/kg K.
rv = 461.5  # gas constant for water vapor. J/kg K.
cp = 1004.  # specific heat at constant pressure (dry air). J/kgK
rlv = 2.53e6  # latent heat for vaporisation
grav = 9.81  # gravity acceleration. m/s^2
mair = 28.967  # Molar mass of air g/mol


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
        # many geometries. points now select only one grid index.
        for g in mask_geoms:
            if isinstance(g, shapely.geometry.Point):
                dists = [haversine.haversine((p[0], p[1]), (g.x, g.y)) for p in points]
                result.append(numpy.argmin(dists))
            else:
                for i in range(len(points)):
                    if g.contains(shapely.geometry.Point(points[i])): result.append(i)
        result = list(set(result))  # remove duplicates
        return result


# Links contents of input directory
def link_dir(files, workdir):
    os.makedirs(workdir)
    for f in files:
        fname = os.path.basename(f)
        os.symlink(os.path.abspath(f), os.path.join(workdir, fname))
