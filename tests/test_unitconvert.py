# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman.unitconvert as uc

class Test_unitconvert:
    def test_build_unit(self):
        assert pytest.approx(uc.unit['mm']) == 10000000.0
        uc.unit['mm'] = 5253.
        assert pytest.approx(uc.unit['mm']) == 5253.
        uc.build_unit()
        assert pytest.approx(uc.unit['mm']) == 10000000.0
        
    def test_set_and_get_in_units(self):
        newton = uc.set_in_units(1e5, 'dyn')
        assert pytest.approx(uc.get_in_units(newton, 'kg*m/s^2')) == 1.0

    def test_set_literal(self):
        value = uc.set_literal('1.124 nm')
        assert pytest.approx(value) == 11.24

    def test_scalar_model(self):
        unit = 'mJ/s^2'
        v = 1234.214
        value = uc.set_in_units(v, unit)
        model = uc.model(value, unit)
        value2 = uc.value_unit(model)
        assert pytest.approx(value) == value2

    def test_vector_model(self):
        unit = 'mJ/s^2'
        v = np.array([1234.214, 346.23])
        value = uc.set_in_units(v, unit)
        model = uc.model(value, unit)
        value2 = uc.value_unit(model)
        assert np.allclose(value, value2)

    def test_tensor_model(self):
        unit = 'mJ/s^2'
        v = np.array([[1234.214, 346.23],[109.124, 235.781]])
        value = uc.set_in_units(v, unit)
        model = uc.model(value, unit)
        value2 = uc.value_unit(model)
        assert np.allclose(value, value2)
