import pytest
import atomman as am

class Test_atomic_mass:
    def test_Fe(self):
        assert am.tools.atomic_mass('Fe') == 55.845
        
    def test_77(self):
        assert am.tools.atomic_mass(77) == 192.217
        
    def test_Not(self):
        with pytest.raises(KeyError):
            mass = am.tools.atomic_mass('Not')
        
    def test_0(self):
        with pytest.raises(KeyError):
            mass = am.tools.atomic_mass(0)
            
    def test_118(self):
        with pytest.raises(ValueError):
            mass = am.tools.atomic_mass(118)