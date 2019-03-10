import numpy
import time
import threading
import netCDF4
import logging
from omuse.units import units

# open a netcdf file for storing lwp fields and vertical profiles
# needs the oifs instance and one les instance for axis information.
# note: currently assumes all les instances have the same parameters
#   - if not, move the axes into the dales subgroups.
#
# for multiple les instances: les output is stored in subgroups.
# if les is None, the les-specific variables and axes are not made.
#

# Logger
log = logging.getLogger(__name__)

# Netcdf root key
cdf_root = None

# Netcdf lock variable
cdf_lock = threading.Lock()  # lock protecting the netCDF root

# Current time step being written
cdf_step = -1

# dict mapping column index to netCDF group handle
# for the extra output columns
output_column_cdf = {}


# Initializes netcdf; when called multiple times, skips initialization
# output_columns is an optional list of extra columns for which output in the netCDF
# is wanted, even though they do not have an embedded LES.
def init_netcdf(nc_name, oifs, les_models, datetime, output_columns=None, append=False, with_surf_vars=True):
    global cdf_root, output_column_cdf
    extra_cols = [] if output_columns is None else output_columns

    if cdf_root:
        cdf_root.close()
        
    if append:
        cdf_root = netCDF4.Dataset(nc_name, "a")
#        print (cdf_root.groups)
        for les in les_models:
            les.cdf = cdf_root.groups[str(les.grid_index)]
        for c in extra_cols:
            idx = c[0]
            cdf = cdf_root.groups[str(idx)]
            print "Warning: untested restart with extra output columns"
            output_column_cdf[idx] = cdf
    else:
        cdf_root = open_netcdf(nc_name, oifs, les_models[0] if any(les_models) else None, datetime)
        for les in les_models:
            les.cdf = create_netcdf_les_subgroup(cdf_root, les, with_surf_vars=with_surf_vars)
        for c in extra_cols:
            idx = c[0]
            lat = c[1]
            lon = c[2]
            cdf = create_netcdf_subgroup(cdf_root, idx, lat, lon, with_surf_vars=with_surf_vars)
            output_column_cdf[idx] = cdf
    return cdf_root


# Updates NetCDF time variable (unlimited) with new time (in s)
def update_time(t):
    global cdf_root, cdf_step
    cdf_step = cdf_root.variables["Time"].shape[0]
    cdf_root.variables["Time"][cdf_step] = t.value_in(units.s)
    log.info("update_time(): step %4d, time %6d s"%(cdf_step, t.value_in(units.s)))


# Flushes the netcdf buffer to disc within thread lock
def sync_root():
    global cdf_root, cdf_lock
    cdf_lock.acquire()
    start = time.time()
    if cdf_root:
        cdf_root.sync()
    cdf_lock.release()
    walltime = time.time() - start
    log.info("netcdf.sync() - %3.1f s" % walltime)


# Opens the new netCDF file
def open_netcdf(nc_name, oifs, les, start_time):
    log.info("opening netcdf %s" % nc_name)
    root_group = netCDF4.Dataset(nc_name, "w")

    # dimensions
    if les:
        root_group.createDimension("x", les.get_itot())
        root_group.createDimension("y", les.get_jtot())
        root_group.createDimension("zf", les.get_ktot())

    root_group.createDimension("oifs_height", oifs.get_ktot())  # height dimension for oifs - note varying pressures and
    # heights
    root_group.createDimension("Time", None)  # time dimension for openIFS steps and forcings

    # coordinate variables
    if les:
        dx = les.dx.value_in(units.m)
        xs = root_group.createVariable("x", "f4", ("x",))
        xs[:] = numpy.linspace(dx / 2, les.xsize.value_in(units.m) - dx / 2, les.get_itot())
        xs.units = 'm'

        dy = les.dy.value_in(units.m)
        ys = root_group.createVariable("y", "f4", ("y",))
        ys[:] = numpy.linspace(dy / 2, les.ysize.value_in(units.m) - dy / 2, les.get_jtot())
        ys.units = 'm'

        zfs = root_group.createVariable("zf", "f4", ("zf",))  # full z levels
        zfs[:] = les.zf.value_in(units.m)
        zfs.units = 'm'

    times = root_group.createVariable("Time", "f4", ("Time",))
    times.units = 's since ' + str(start_time)

    return root_group


# Creates subgroup in the netcdf file for the given les instance
# add variables to the subgroup that are LES-specific.
# Variables from openIFS are added in create_netcdf_subgroup()
def create_netcdf_les_subgroup(rootgrp, les, with_surf_vars=True):
    subgroup = str(les.grid_index)
    grp = create_netcdf_subgroup(rootgrp, subgroup, les.lat, les.lon, with_surf_vars)

    # vertical profiles - using absolute heights from Dales, zf height dimension
    # using the Time dimension - store once every large step
    for name, unit in (('u', 'm/s'),  # obtained from dales
                       ('v', 'm/s'),
                       ('thl', 'K'),
                       ('qt', '1'),
                       ('ql', '1'),
                       ('ql_ice', '1'),
                       ('ql_water', '1'),
                       ('qr', '1'),
                       ('t', 'K'),  # temperature calculated in spifs using OpenIFS pressures
                       ('t_', 'K'),  # temperature calculated in dales
                       ('f_u', 'm/s'),  # forcings on Dales
                       ('f_v', 'm/s'),
                       ('f_thl', 'K/s'),
                       ('f_qt', '1/s'),
                       ('presf', 'Pa/s'),
                       ('qt_std', '1'),
                       ('qt_alpha', '1/s'),
                       ('qt_beta', '1')):
        p = grp.createVariable(name, "f4", ("Time", "zf"))
        p.units = unit

    # vertical profiles - using heights from OpenIFS
    for name, unit in (
            ('f_U', 'm/s'),  # forcings on OpenIFS
            ('f_V', 'm/s'),
            ('f_T', 'K/s'),
            ('f_SH', '1/s'),
            ('f_QL', '1/s'),
            ('f_QI', '1/s'),
            ('f_A', '1/s')):
        p = grp.createVariable(name, "f4", ("Time", "oifs_height"))
        p.units = unit

    return grp


# Creates a subgroup in the netcdf file
def create_netcdf_subgroup(rootgrp, subgroup, lat, lon, with_surf_vars=True):
    grp = rootgrp.createGroup(str(subgroup))

    # vertical profiles - using heights from OpenIFS
    for name, unit in (('U', 'm/s'),
                       ('V', 'm/s'),
                       ('T', 'K'),
                       ('SH', '1'),
                       ('QL', '1'),
                       ('QI', '1'),
                       ('Pf', 'Pa'),
                       ('Ph', 'Pa'),
                       ('Tv', 'K'),
                       ('Zf', 'm'),
                       ('Zh', 'm'),
                       ('THL', 'K'),
                       ('QT', '1'),
                       ('A', '1')):
        p = grp.createVariable(name, "f4", ("Time", "oifs_height"))
        p.units = unit

    # Surface fields (LES scalars)
    srf = [('Psurf',     'Pa'),
           ('rain',      'kg / m^2'),
           ('rainrate',  'kg / m^2h')
          ]

    if with_surf_vars:
        srf += [('z0m', 'm')]
        srf += [('z0h', 'm')]
        srf += [('wthl', 'K m/s')]
        srf += [('wqt', 'kg m/s')]
        srf += [('TLflux', 'W/m^2')]
        srf += [('TSflux', 'W/m^2')]
        srf += [('SHflux', 'kg / m^2s')]
        srf += [('QLflux', 'kg / m^2s')]
        srf += [('QIflux', 'kg / m^2s')]

    for name, unit in srf:
        p = grp.createVariable(name, "f4", ("Time",))
        p.units = unit

    # coordinates of this les instance.
    # NOTE for radiation etc in Dales, need to set these coordinates
    # for the worker code using the interface - this is not yet done
    lat_ = grp.createVariable("lat", "f4")
    lat_.units = 'deg'
    lon_ = grp.createVariable("lon", "f4")
    lon_.units = 'deg'
    lat_[:] = lat
    lon_[:] = lon
    return grp


def write_les_data(les, **kwargs):
    global cdf_lock
    lock = kwargs.get("lock", False)
    if lock:
        cdf_lock.acquire()
    for var, arr in kwargs.iteritems():
        if var == "lock":
            continue  # variable argument list nonsense
        ncvar = les.cdf.variables.get(var, None)
        if ncvar:
            ncvar[cdf_step] = kwargs[var]
        else:
            log.error("Attempt to write profile to uninitialized variable %s" % var)
    if lock:
        cdf_lock.release()


# alternative version used for the extra output columns
# todo: merge with write_les_data ?
# can also store the les columns col_index->cdf mapping in the dict
def write_netCDF_data(column_index, **kwargs):
    global cdf_lock

    cdf_handle = output_column_cdf[column_index]
    lock = kwargs.get("lock", False)
    if lock:
        cdf_lock.acquire()
    for var, arr in kwargs.iteritems():
        if var == "lock": continue  # variable argument list nonsense
        ncvar = cdf_handle.variables.get(var, None)

        if ncvar:
            try:
                ncvar[cdf_step] = arr
            except IndexError:
                log.error("write_netCDF_data: Variable %s has length %d, should be %d" % (var, len(arr), len(ncvar[cdf_step])))
            except ValueError:
                log.error("write_netCDF_data: Variable %s has a weird type"%var)
                log.error(str(var.unit)) # try to write the unit
        else:
            log.error("Attempt to write profile to uninitialized variable %s" % var)
    if lock:
        cdf_lock.release()
