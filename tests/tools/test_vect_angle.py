# coding: utf-8

# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# https://docs.pytest.org/en/latest/
import pytest

from atomman.tools import vect_angle

class Test_vect_angle:
    def test_0(self):
        assert pytest.approx(vect_angle([1,0,0], [1,0,0]), 0.0)
        
    def test_180(self):
        assert pytest.approx(vect_angle([1,0,0], [-1,0,0]), 180.0)
        
    def test_90(self):
        assert pytest.approx(vect_angle([1,0,0], [0,1,0]), 90.0)
        
    def test_45(self):
        assert pytest.approx(vect_angle([1,1,0], [1,1,0]), 45.0)