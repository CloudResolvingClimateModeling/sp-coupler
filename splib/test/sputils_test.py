import numpy
from splib import sputils
from splib import spdummy

class Testsputils(object):

    tolerance = 1.e-10


    def test_rms(self):
        a = 3.204
        b = -1.2092
        c = 9.6231
        numbers = [a,b,c]
        assert abs(sputils.rms(numpy.array(numbers)) - numpy.sqrt((a * a + b * b + c * c)/3)) < self.tolerance


    def test_rms_repeat(self):
        a = 3.204
        n = 23
        numbers = [a for i in range(n)]
        assert abs(sputils.rms(numpy.array(numbers)) - a) < self.tolerance


    def test_exner(self):
        a = 2.03947
        p = a * sputils.pref0
        assert abs(numpy.log(sputils.exner(p)) - numpy.log(a) * sputils.rd / sputils.cp) < self.tolerance


    def test_exner_unity(self):
        p = sputils.pref0
        assert abs(sputils.exner(p) - 1) < self.tolerance


    def test_iexner(self):
        a = 12.03947
        p = a * sputils.pref0
        assert abs(sputils.exner(p)*sputils.iexner(p) - 1) < self.tolerance


    def test_get_closest_points(self):
        points = [(52.314970,4.824198),(52.379932,4.897997),(52.387264,5.082968),(52.278097,5.021635)]
        target = (52.356591, 4.954541)
        assert sputils.find_closest_points(points,target)[0] == 1
