# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am
import atomman.unitconvert as uc

class Test_ElasticConstants:
    def test_default(self):
        c = am.ElasticConstants()
        assert np.allclose(c.Cij, np.zeros((6,6)))
        assert np.allclose(c.Cij9, np.zeros((9,9)))
        assert np.allclose(c.Cijkl, np.zeros((3,3,3,3)))