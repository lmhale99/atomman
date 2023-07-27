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

    def test_isotropic(self):
        pass

    def test_cubic(self):

        C11 = uc.set_in_units(140, 'GPa')
        C12 = uc.set_in_units(100, 'GPa')
        C44 = uc.set_in_units(90, 'GPa')

        C = am.ElasticConstants(C11=C11, C12=C12, C44=C44)
        assert C.is_normal('cubic')

        Cij = np.array([
            [C11, C12, C12, 0.0, 0.0, 0.0],
            [C12, C11, C12, 0.0, 0.0, 0.0],
            [C12, C12, C11, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, C44, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, C44, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, C44],
        ])
        assert np.allclose(C.Cij, Cij)
