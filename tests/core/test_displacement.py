# coding: utf-8

# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from copy import deepcopy

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am

def test_displacement():
    a = 3.18
    c = 5.17
    ucell = am.System(atoms=am.Atoms(pos=[[0.0, 0.0, 0.0],
                                        [1/3, 2/3, 0.5]]), 
                    box=am.Box.hexagonal(a, c),
                    scale=True)
    system_0 = ucell.supersize(10, 10, 10)
    system_1 = deepcopy(system_0)

    system_1.atoms.pos[:] += 5
    system_1.wrap()

    assert not np.allclose(system_1.atoms.pos - system_0.atoms.pos, 5.0)
    assert np.allclose(am.displacement(system_0, system_1), 5.0)