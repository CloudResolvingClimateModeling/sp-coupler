from __future__ import division
import logging
import numpy
import datetime
from amuse.community import units

# Logger
log = logging.getLogger(__name__)


# Base class for dummy models
# noinspection PyMethodMayBeStatic
class dummy_base(object):

    def __init__(self, nprocs):
        self.step = 0
        self.timestep = 600 | units.s
        self.starttime = datetime.datetime(2000, 1, 1)
        self.model_time = 0 | units.s
        self.number_of_workers = nprocs
        self.support_async = False
        self.itot, self.jtot, self.ktot = 10, 10, 20

    def get_itot(self):
        return self.itot

    def get_jtot(self):
        return self.jtot

    def get_ktot(self):
        return self.ktot

    def get_timestep(self):
        return self.timestep

    def get_model_time(self):
        return self.model_time

    def evolve_model_single_step(self):
        self.model_time += self.timestep
        return True

    def evolve_model_until_cloud_scheme(self):
        return True

    def evolve_model_cloud_scheme(self):
        return True

    def evolve_model_from_cloud_scheme(self):
        self.model_time += self.timestep
        return True

    def cleanup_code(self):
        return True

    def stop(self):
        return True


# Dummy climate model class
# noinspection PyMethodMayBeStatic
class dummy_gcm(dummy_base):

    def __init__(self, nprocs):
        super(dummy_gcm, self).__init__(nprocs)
        self.time_per_gridpoint = 0.0001
        self.num_lats = 20
        self.latitudes = numpy.empty([self.num_lats], dtype=numpy.float64)
        self.num_lons = 40
        self.longitudes = numpy.empty([self.num_lons], dtype=numpy.float64)
        self.ktot = 20
        self.mask = set([])
        self.step_time = 0.
        log.info("Initialized dummy gcm with %d latitudes, %d longitudes and %d vertical layers" % (
            self.num_lats, self.num_lons, self.ktot))

    def get_itot(self):
        return self.num_lons

    def get_jtot(self):
        return self.num_lats

    def get_ktot(self):
        return self.ktot

    def get_start_datetime(self):
        return self.starttime

    def initialize_code(self):
        log.info("Initialized dummy gcm code")

    def commit_parameters(self):
        self.step_time = self.timestep * int(self.num_lats * self.num_lons * self.time_per_gridpoint)
        log.info("Initialized dummy gcm parameters")

    def commit_grid(self):
        lats = numpy.fromfunction(lambda i: 180 * (i / self.num_lats) - 90., (self.num_lats,))
        lons = numpy.fromfunction(lambda i: 360 * (i / self.num_lons), (self.num_lons,))
        self.latitudes = numpy.repeat(lats, len(lons)) | units.deg
        self.longitudes = numpy.tile(lons, len(lats)) | units.deg
        log.info("Initialized dummy gcm grid")

    def set_mask(self, i):
        self.mask.add(i)

    def get_field(self, name, i, k):
        func_hor, func_vert, norm, unit = self.field_helper(name)
        field = norm * func_hor(self.latitudes[i], self.longitudes[i]) * func_vert(k)
        return field | unit if unit else field

    def get_profile_field(self, name, index):
        func_hor, func_vert, norm, unit = self.field_helper(name)
        func = numpy.vectorize(lambda i: norm * func_hor(self.latitudes[i], self.longitudes[i]))
        factors = func(index)
        values = self.ktot + 1 if name in ["Phalf", "Zh"] else self.ktot
        field = numpy.outer(factors, numpy.fromfunction(func_vert, (values,)))
        return field | unit if unit else field

    def get_profile_fields(self, name, index):
        func_hor, func_vert, norm, unit = self.field_helper(name)
        func = numpy.vectorize(lambda i: norm * func_hor(self.latitudes[i], self.longitudes[i]))
        factors = func(index)
        values = self.ktot + 1 if name in ["Phalf", "Zh"] else self.ktot
        field = numpy.outer(factors, numpy.fromfunction(func_vert, (values,)))
        return field | unit if unit else field

    def get_volume_field(self, name):
        func_hor, func_vert, norm, unit = self.field_helper(name)
        values = self.ktot + 1 if name in ["Phalf", "Zh"] else self.ktot
        field = norm * numpy.fromfunction(lambda i, k: func_hor(self.latitudes[i], self.longitudes[i]) * func_vert(k),
                                          (self.num_lats * self.num_lons, values,))
        return field | unit if unit else field

    def get_layer_field(self, name, index):
        func_hor, func_vert, norm, unit = self.field_helper(name)
        factor = norm * func_vert(index)
        field = factor * numpy.fromfunction(lambda i: func_hor(self.latitudes[i], self.longitudes[i]),
                                            (self.num_lats * self.num_lons,))
        return field | unit if unit else field

    def field_helper(self, name):
        func_xy = lambda x, y: 1 + 0.3 * numpy.cos(x * y)
        func_k = lambda k: 1. / ((k - 0.5 * self.ktot) ** 2 + 1)
        norm = 1
        unit = None
        if name in ['U', 'V', 'W']:
            norm = 10.
            unit = units.m / units.s
        elif name in ['T']:
            norm = 273.
            unit = units.K
        elif name in ["Pfull"]:
            norm = 100000.
            func_xy = lambda x, y: 1.
            func_k = lambda k: numpy.exp(-4 * (self.ktot - k - 0.5) / self.ktot)
            unit = units.Pa
        elif name in ["Phalf"]:
            norm = 100000.
            func_xy = lambda x, y: 1.
            func_k = lambda k: numpy.exp(-4 * (self.ktot - k) / self.ktot)
            unit = units.Pa
        elif name in ["Zf"]:
            norm = 50000.
            func_xy = lambda x, y: 1.
            func_k = lambda k: (self.ktot - k - 0.5) / self.ktot
            unit = units.m
        elif name in ["Zh"]:
            norm = 50000.
            func_xy = lambda x, y: 1.
            func_k = lambda k: (self.ktot - k) / self.ktot
            unit = units.m
        return func_xy, func_k, norm, unit

    def set_profile_tendency(self, field, index, vals):
        log.info("Setting profile tendency for %s at grid point %d" % (field, index))

    def set_vdf_in_sp_mask(self, value):
        log.info("Setting vdf process switch to %s" % str(value))


# Dummy LES model class
# noinspection PyUnresolvedReferences,PyMethodMayBeStatic,PyPep8Naming
class dummy_les(dummy_base):

    def __init__(self, nprocs):
        super(dummy_les, self).__init__(nprocs)
        self.time_per_gridpoint = 0.0001
        self.itot = 8
        self.jtot = 8
        self.k = 20
        self.dx = 100 | units.m
        self.dy = 100 | units.m
        self.dz = 200 | units.m
        self.xsize = self.itot * self.dx
        self.ysize = self.jtot * self.dy
        self.zsize = self.k * self.dz
        self.zf = numpy.empty([self.k])
        self.zh = numpy.empty([self.k])
        self.sp = 100000. | units.Pa
        self.step_time = 0.
        log.info("Initialized dummy les with %d x-coords, %d y-coords and %d z-coords" % (self.itot, self.jtot, self.k))

    def get_itot(self):
        return self.itot

    def get_jtot(self):
        return self.jtot

    def get_ktot(self):
        return self.k

    # noinspection PyMethodMayBeStatic
    def initialize_code(self):
        log.info("Initialized dummy les code")

    def commit_parameters(self):
        self.step_time = self.itot * self.jtot * self.k * self.time_per_gridpoint
        log.info("Initialized dummy les parameters")

    def commit_grid(self):
        self.zf = numpy.fromfunction(lambda k: k * self.dz, (self.k,), dtype=int)
        self.zh = numpy.fromfunction(lambda k: (k + 0.5) * self.dz, (self.k,), dtype=int)
        log.info("Initialized dummy les grid")

    def evolve_model(self, stop_time, exactEnd):
        if exactEnd:
            log.info("Evolving LES model exactly to time %s" % str(stop_time))
        else:
            log.info("Evolving LES model to at least time %s" % str(stop_time))
        self.model_time = stop_time

    # noinspection PyMethodMayBeStatic
    def set_field(self, fid, values):
        log.info("Setting vertical profile for field %s with values %s" % (fid, str(values)))

    def set_surface_pressure(self, value):
        log.info("Setting surface pressure to %f Pa" % value.value_in(units.Pa))
        self.sp = value

    def get_surface_pressure(self):
        return self.sp

    def get_field(self, name):
        if name == "TWP":
            return numpy.fromfunction(lambda i, j: numpy.sin(6.28 / (i + j + 1)) + numpy.cos(6.28 / (i + j + 1)),
                                      (self.itot, self.jtot))
        if name == "LWP":
            return numpy.fromfunction(lambda i, j: numpy.sin(6.28 / (i + j + 1)), (self.itot, self.jtot))
        if name == "RWP":
            return numpy.fromfunction(lambda i, j: numpy.cos(6.28 / (i + j + 1)), (self.itot, self.jtot))
        return None

    def get_profile_field(self, name):
        log.info("Getting LES profile for variable %s" % name)
        if name in ["U", "V", "W"]:
            return numpy.sin(
                6.28 * self.zf.value_in(units.m) / (self.dz.value_in(units.m) * self.k)) | units.m / units.s
        if name in ["THL", "T"]:
            return (283.0 + 10. * numpy.cos(
                6. * self.zf.value_in(units.m) / (self.dz.value_in(units.m) * self.k))) | units.K
        if name in ["QT", "A"]:
            return 0.5 + 0.2 * numpy.cos(6. * self.zf.value_in(units.m) / (self.dz.value_in(units.m) * self.k))
        if name in ["QL"]:
            return 0.5 + 0.2 * numpy.sin(6. * self.zf.value_in(units.m) / (self.dz.value_in(units.m) * self.k))
        if name in ["QR"]:
            return 0.0001 * numpy.sin(6. * self.zf.value_in(units.m) / (self.dz.value_in(units.m) * self.k))
        if name in ["zh"]:
            return self.zh
        if name in ["zf"]:
            return self.zf
        if name in ["ph"]:
            return 100000. * numpy.exp(-self.zh.value_in(units.m)) | units.Pa
        if name in ["pf"]:
            return 100000. * numpy.exp(-self.zf.value_in(units.m)) | units.Pa
        return None

    def get_profile_U(self):
        return self.get_profile_field("U")

    def get_profile_V(self):
        return self.get_profile_field("V")

    def get_profile_W(self):
        return self.get_profile_field("W")

    def get_profile_T(self):
        return self.get_profile_field("T")

    def get_profile_QT(self):
        return self.get_profile_field("QT")

    def get_profile_QL(self):
        return self.get_profile_field("QL")

    def get_profile_QL_ice(self):
        return 0.1 * self.get_profile_field("QL")

    def get_profile_QL_water(self):
        return 0.9 * self.get_profile_field("QL")

    def get_profile_QR(self):
        return self.get_profile_field("QR")

    def get_profile_THL(self):
        return self.get_profile_field("THL")

    def get_presf(self):
        return self.get_profile_field("pf")

    def get_zf(self):
        return self.get_profile_field("zf")

    def get_presh(self):
        return self.get_profile_field("ph")

    def get_zh(self):
        return self.get_profile_field("zh")

    def get_cloudfraction(self, i):
        indices = numpy.clip(i, 0, self.k - 1)
        return self.get_profile_field("A")[indices]

    def set_tendency_U(self, values):
        log.info("Setting U-tendency to %s" % str(values))

    def set_tendency_V(self, values):
        log.info("Setting V-tendency to %s" % str(values))

    def set_tendency_THL(self, values):
        log.info("Setting THL-tendency to %s" % str(values))

    def set_tendency_QT(self, values):
        log.info("Setting QT-tendency to %s" % str(values))

    def set_tendency_surface_pressure(self, values):
        log.info("Setting SP-tendency to %s" % str(values))

    def set_multiplicative_qt_forcing(self, values):
        log.info("Setting multiplicative qt forcing to %s" % str(values))

    def set_fluctuation_forcing(self, values):
        log.info("Setting variance qt forcing to %s" % str(values))

    def set_ref_profile_QL(self, values):
        log.info("Setting reference ql profile to %s" % str(values))
