import pytest
import atomman.unitconvert as uc

class Test_unitconvert:
    def test_unit_get(self):
        assert pytest.approx(uc.unit['mm'], 10000000.0)

    def test_unit_set(self):
        uc.unit['mm'] = 5253.
        assert pytest.approx(uc.unit['mm'], 5253.)

    def test_build_unit(self):
        uc.unit['mm'] = 5253.
        uc.build_unit()
        assert pytest.approx(uc.unit['mm'], 10000000.0)
        
    def test_newton(self):
        newton = uc.set_in_units(1e5, 'dyn')
        assert pytest.approx(uc.get_in_units(newton, 'kg*m/s^2'), 1.0)