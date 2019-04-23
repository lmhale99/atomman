# coding: utf-8

# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# https://docs.pytest.org/en/latest/
import pytest

from atomman.tools import atomic_mass, atomic_symbol, atomic_number

class Test_atomic_info:
    def test_atomic_symbol(self):
        assert atomic_symbol(74) == 'W'
        
        with pytest.raises(IndexError):
            atomic_symbol(0)
        
        with pytest.raises(IndexError):
            atomic_symbol(119)

    def test_atomic_number(self):
        assert atomic_number('Fe') == 26
        
        with pytest.raises(ValueError):
            atomic_number('Nt')
        
    def test_atomic_mass(self):
        assert atomic_mass(77) == 192.217

        assert atomic_mass('He-9') == 9.043946

        with pytest.raises(ValueError):
            atomic_mass(118)

        with pytest.raises(IndexError):
            atomic_mass(119)