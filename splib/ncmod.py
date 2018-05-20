from __future__ import division

import logging

import netCDF4
import numpy
import parser
from amuse.community import units
from amuse.community.interface.common import CommonCode

# Logger
log = logging.getLogger(__name__)


# Base class for netcdf-based models
# noinspection PyMethodMayBeStatic
class netcdf_model_base(CommonCode):

    def __init__(self, ncfile):
        self.dataset = netCDF4.Dataset(ncfile)
        self.model_time, self.time_step = self.get_time_info(self.dataset)
        self.step = 0
        self.number_of_workers = 1
        self.support_async = False

    def __del__(self):
        self.dataset.close()

    # Returns start time and time step from nc dataset
    @staticmethod
    def get_time_info(ds):
        tim_var = ds.variables["Time"]
        tim_vals = tim_var[:]
        if len(tim_vals) == 0:
            raise Exception("Zero time values found in spifs.nc file")
        if len(tim_vals) == 1:
            return 0 | units.s, tim_vals[0] | units.s
        dt = (tim_vals[1] - tim_vals[0]) | units.s
        return tim_vals[0] | units.s - dt, dt

    @staticmethod
    def get_units(var):
        units = getattr(var, "units", None)
        if not units or units == "1":
            return None
        # noinspection PyUnusedLocal
        m, s, K, Pa, W = units.m, units.s, units.K, units.Pa, units.W
        expression = parser.expr(units).compile()
        result = eval(expression)
        return result

    def initialize_code(self):
        log.error("Initializing code netcdf in abstract base class")

    def commit_parameters(self):
        log.error("Committing paramters  netcdf in abstract base class")

    def get_timestep(self):
        return self.time_step

    def get_model_time(self):
        return self.model_time

    def evolve_model_single_step(self):
        self.step += 1
        self.model_time += self.time_step
        return True

    def evolve_model_until_cloud_scheme(self):
        return True

    def evolve_model_cloud_scheme(self):
        return True

    def evolve_model_from_cloud_scheme(self):
        self.step += 1
        self.model_time += self.time_step
        return True

    def cleanup_code(self):
        self.dataset.close()
        return True

    def stop(self):
        return True


# GCM based on the spifs.nc output file. It has its grid equal to the
# superparametrization mask
class netcdf_gcm(netcdf_model_base):

    def __init__(self, ncfile):
        super(netcdf_gcm, self).__init__(ncfile)
        self.latitudes, self.longitudes, self.group_names = [], [], []
        pts = len(self.dataset.groups.keys())
        self.num_lats = pts
        self.num_lons = pts
        self.ktot = len(self.dataset.dimensions["oifs_height"][:])
        self.mask = set(range(pts))
        log.info("Initialized netcdf gcm with %d latitudes, %d longitudes and %d vertical layers" % (
            self.num_lats, self.num_lons, self.ktot))

    def initialize_code(self):
        log.info("Initialize netcdf gcm code")

    def commit_parameters(self):
        log.info("Commit netcdf gcm parameters")

    def commit_grid(self):
        for k, v in self.dataset.groups:
            self.group_names.append(k)
            self.latitudes.append(v.variables["latitude"][0])
            self.longitudes.append(v.variables["longitude"][0])
        log.info("Initialized netcdf gcm grid")

    def set_mask(self, i):
        if i < self.num_lats:
            self.mask.add(i)
        else:
            log.error("Attempt to add index %d to sp mask, which is out of range for this netcdf file-based model" % i)

    def get_field(self, name, i, k):
        group_name = self.group_names[i]
        vardict = self.dataset.groups[group_name].variables[name]
        values = vardict[self.step, k]
        units = netcdf_model_base.get_units(vardict)
        return values | units if units else values

    def get_profile_field(self, name, index):
        group_name = self.group_names[index]
        vardict = self.dataset.groups[group_name].variables[name]
        values = vardict[self.step, :]
        units = netcdf_model_base.get_units(vardict)
        return values | units if units else values

    def get_profile_fields(self, name, indices):
        group_names = self.group_names[indices]
        if not any(group_names):
            return []
        values = numpy.vstack([self.dataset.groups[g].variables[name][self.step, :] for g in group_names])
        units = netcdf_model_base.get_units(self.dataset.groups[group_names[0]].variables[name])
        return values | units if units else values

    def get_volume_field(self, name):
        if not any(self.group_names):
            return []
        values = numpy.vstack([self.dataset.groups[g].variables[name][self.step, :] for g in self.group_names])
        units = netcdf_model_base.get_units(self.dataset.groups[self.group_names[0]].variables[name])
        return values | units if units else values

    def get_layer_field(self, name, index):
        if not any(self.group_names):
            return []
        values = numpy.array([self.dataset.groups[g].variables[name][self.step, index] for g in self.group_names])
        units = netcdf_model_base.get_units(self.dataset.groups[self.group_names[0]].variables[name])
        return values | units if units else values

    def set_profile_tendency(self, field, index, vals):
        log.info("Setting profile tendency for %s at grid point %d" % (field, index))
        ncprof = self.get_profile_field(field, index)
        log.info("Difference with values in file are %s" % str(vals - ncprof))


# LES based on the spifs.nc output file
# noinspection PyPep8Naming
class netcdf_les(netcdf_model_base):

    def __init__(self, ncfile, group_index=-1):
        super(netcdf_les, self).__init__(ncfile)
        self.itot = self.dataset.dimensions["x"].shape[0]
        self.jtot = self.dataset.dimensions["y"].shape[0]
        self.zf = numpy.copy(self.dataset.variables["zf"][:])
        self.k = len(self.zf)
        self.xsize = self.dataset.variables["x"][-1] - self.dataset.variables["x"][0]
        self.ysize = self.dataset.variables["y"][-1] - self.dataset.variables["y"][0]
        self.zsize = self.zf[-1] - self.dataset.variables["zf"][0]
        self.zh = 0.5*(self.zf[1:] + self.zf[:-1])
        self.group_index = group_index
        if group_index >= 0:
            self.group = self.dataset.groups[group_index]
            self.sp = self.group.variables["Psurf"][0] | units.Pa
        self.step_time = 0.
        log.info("Initialized netcdf les with %d x-coords, %d y-coords and %d z-coords" % (self.itot, self.jtot, self.k))

    @property
    def grid_index(self):
        return self.group_index

    @grid_index.setter
    def grid_index(self, i):
        self.group_index = i
        if self.group_index > i:
            self.group = self.dataset.groups[self.group_index]
            self.sp = self.group.variables["Psurf"][0] | units.Pa

    def initialize_code(self):
        log.info("Initialize netcdf les code")

    def commit_parameters(self):
        log.info("Commit netcdf les parameters")

    def commit_grid(self):
        lat = getattr(self, "lat", None)
        nclat = self.group.variables["lat"][0]
        if lat is None:
            setattr(self, "lat", nclat)
        if lat != nclat:
            log.warning("Latitude mismatch: this model has been assigned a latitude %f whereas the netcdf file "
                        "locates it at latitude %f" % (lat, nclat))
        lon = getattr(self, "lon", None)
        nclon = self.group.variables["lon"][0]
        if lon is None:
            setattr(self, "lon", nclon)
        if lon != nclon:
            log.warning("Longitude mismatch: this model has been assigned a longitude %f whereas the netcdf file "
                        "locates it at longitude %f" % (lon, nclon))
        log.info("Initialized netcdf les grid")

    def evolve_model(self, stop_time, exact_end):
        if exact_end:
            log.info("Evolving les model exactly to time %s" % str(stop_time))
        else:
            log.info("Evolving les model to at least time %s" % str(stop_time))
        self.model_time = stop_time

    # noinspection PyMethodMayBeStatic
    def set_field(self, fid, values):
        log.info("Setting vertical profile for field %s with values %s" % (fid, str(values)))
        nc_values = self.group.variables[fid][self.step,:]
        if values != nc_values:
            log.warning("The difference with the netcdf file is %s" % str(abs(values - nc_values)))

    def set_surface_pressure(self, value):
        log.info("Setting surface pressure to %f Pa" % value.value_in(units.Pa))
        nc_sp = self.group.variables["Psurf"][self.step]
        if value != nc_sp:
            log.warning("The difference with the netcdf file is %d" % abs(value - nc_sp))

    def set_tendency(self, name, values):
        log.info("Setting %s-tendency to %s" % (name, str(values)))
        nc_values = self.group.variables["f_" + "name"][self.step, :]
        if values != nc_values:
            log.warning("The difference with the netcdf file is %d" % abs(values - nc_values))

    def get_surface_pressure(self):
        return self.group.variables["Psurf"][self.step]

    def get_field(self, name):
        return self.group.variables[name][self.step]

    def get_profile_field(self, name):
        return self.group.variables[name][self.step, :]

    def get_profile_U(self):
        return self.get_profile_field("u")

    def get_profile_V(self):
        return self.get_profile_field("v")

    def get_profile_W(self):
        return self.get_profile_field("w")

    def get_profile_T(self):
        return self.get_profile_field("t_")

    def get_profile_QT(self):
        return self.get_profile_field("qt")

    def get_profile_QL(self):
        return self.get_profile_field("ql")

    def get_profile_QL_ice(self):
        return self.get_profile_field("ql_ice")

    def get_profile_QL_water(self):
        return self.get_profile_field("ql_water")

    def get_profile_QR(self):
        return self.get_profile_field("qr")

    def get_profile_THL(self):
        return self.get_profile_field("thl")

    def get_cloudfraction(self, indices):
        A = self.get_profile_field("A")
        if len(indices) == len(A):
            return A[::-1]
        else:
            return numpy.zeros(len(indices))

    def set_tendency_U(self, values):
        self.set_tendency("u", values)

    def set_tendency_V(self, values):
        self.set_tendency("v", values)

    def set_tendency_THL(self, values):
        self.set_tendency("thl", values)

    def set_tendency_QT(self, values):
        self.set_tendency("qt", values)

    def set_tendency_surface_pressure(self, values):
        log.info("Setting presf-tendency to %s" % str(values))
        nc_values = self.group.variables["presf"][self.step, :]
        if values != nc_values:
            log.warning("The difference with the netcdf file is %d" % abs(values - nc_values))
