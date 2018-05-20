import numpy
from splib import spcpl
from splib import spdummy

class Testspcpl(object):

    tolerance = 1.e-10


    def test_get_cloudfraction(self):
        les = spdummy.dummy_les(1)
        les.commit_grid()
        les.gcm_Zh = numpy.array([100000.,1000.,100.,10.,1.,0.])
        A = spcpl.get_cloud_fraction(les)
        assert abs(A[0] - (0.5 + 0.2*numpy.cos(6.*(1. - les.k)/les.k))) < self.tolerance
        assert abs(A[-1] - (0.5 + 0.2)) < self.tolerance
