# coding: utf-8

# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am

def test_fcc():
    
    # Build fcc test system
    a = 4.05
    ucell = am.System(atoms=am.Atoms(pos=[[0.0, 0.0, 0.0],
                                          [0.5, 0.5, 0.0],
                                          [0.5, 0.0, 0.5],
                                          [0.0, 0.5, 0.5]]), 
                    box=am.Box.cubic(a),
                    scale=True)
    system = ucell.supersize(20, 20, 20)
    
    cutoff = 0.9 * a

    # Test all periodic boundaries
    system.pbc = [True, True, True]
    neighbors = system.neighborlist(cutoff=cutoff)
    assert np.isclose(neighbors.coord.mean(), 12.0)

    # Test no periodic boundaries
    system.pbc = [False, False, False]
    neighbors = system.neighborlist(cutoff=cutoff)
    assert np.isclose(neighbors.coord.mean(), 11.4075)

def test_bcc():
    
    # Build bcc test system
    a = 2.865
    ucell = am.System(atoms=am.Atoms(pos=[[0.0, 0.0, 0.0],
                                          [0.5, 0.5, 0.5]]), 
                      box=am.Box.cubic(a),
                      scale=True)
    system = ucell.supersize(30, 30, 30)
    
    cutoff = 0.9 * a

    # Test all periodic boundaries
    system.pbc = [True, True, True]
    neighbors = system.neighborlist(cutoff=cutoff)
    assert np.isclose(neighbors.coord.mean(), 8.0)

    # Test no periodic boundaries
    system.pbc = [False, False, False]
    neighbors = system.neighborlist(cutoff=cutoff)
    assert np.isclose(neighbors.coord.mean(), 7.6066)

def test_hcp():
    
    # Build hcp test system
    a = 3.18
    c = 5.17
    ucell = am.System(atoms=am.Atoms(pos=[[0.0, 0.0, 0.0],
                                          [1/3, 2/3, 0.5]]), 
                      box=am.Box.hexagonal(a, c),
                      scale=True)
    system = ucell.supersize(30, 30, 30)
    
    cutoff = 1.1 * a

    # Test all periodic boundaries
    system.pbc = [True, True, True]
    neighbors = system.neighborlist(cutoff=cutoff)
    assert np.isclose(neighbors.coord.mean(), 12.0)

    # Test no periodic boundaries
    system.pbc = [False, False, False]
    neighbors = system.neighborlist(cutoff=cutoff)
    assert np.isclose(neighbors.coord.mean(), 11.4411)