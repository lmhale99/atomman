import pytest
import atomman as am

class Test_atomic_info:
    def test_atomic_symbol(self):
        assert am.tools.atomic_symbol(74) == 'W'
        
        with pytest.raises(IndexError):
            am.tools.atomic_symbol(0)
        
        with pytest.raises(IndexError):
            am.tools.atomic_symbol(119)

    def test_atomic_number(self):
        assert am.tools.atomic_number('Fe') == 26
        
        with pytest.raises(ValueError):
            am.tools.atomic_number('Nt')
        
    def test_atomic_mass(self):
        assert am.tools.atomic_mass(77) == 192.217

        assert am.tools.atomic_mass('He-9') == 9.043946

        with pytest.raises(ValueError):
            am.tools.atomic_mass(118)

        with pytest.raises(IndexError):
            am.tools.atomic_mass(119)