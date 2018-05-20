import numpy
import os
import shapely.geometry
import netCDF4
from amuse.community import units
from splib import splib
from splib import spdummy

class Testsplib(object):

    def test_lwp_after_run(self):
        lesdt,steps = 10,5
        config = {"gcm_type":"dummy","les_type":"dummy","les_dt":lesdt}
        splib.initialize(config,[shapely.geometry.Point(50.0,2.0)])
        splib.run(steps)
        splib.finalize()
        filename = os.path.join(splib.output_dir,splib.output_name)
        output = netCDF4.Dataset(filename,'r')
        for group in output.groups:
            print output.variables.keys(),output[str(group) + "/lwp"].shape
            assert output[str(group) + "/lwp"].shape[0] == steps * splib.gcm_model.get_timestep().value_in(units.s) / lesdt
