# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am
from atomman.tools.miller import (vector3to4, vector4to3,
                                  vector_conventional_to_primitive,
                                  vector_crystal_to_cartesian,
                                  vector_primitive_to_conventional)

def test_vector3to4():
    assert np.all(vector3to4([ 1, 2, 3]) == np.array([ 0, 1,-1, 3]))
    assert np.all(vector3to4([ 0, 3, 0]) == np.array([-1, 2,-1, 0]))
    assert np.all(vector3to4([-3, 0, 1]) == np.array([-2, 1, 1, 1]))
    assert np.all(vector3to4([-2, 2, 2]) == np.array([-2, 2, 0, 2]))
    with pytest.raises(ValueError):
        vector3to4([ 1, 2, 3, 4])


def test_vector4to3():
    assert np.all(vector4to3([ 0, 1,-1, 3]) == np.array([ 1, 2, 3]))
    assert np.all(vector4to3([-1, 2,-1, 0]) == np.array([ 0, 3, 0]))
    assert np.all(vector4to3([-2, 1, 1, 1]) == np.array([-3, 0, 1]))
    assert np.all(vector4to3([-2, 2, 0, 2]) == np.array([-2, 2, 2]))
    with pytest.raises(ValueError):
        vector4to3([ 1, 2, 3])

def test_vector_crystal_to_cartesian():
    # Test 3 indices
    box = am.Box(a=3, b=4, c=5, alpha=90, beta=90, gamma=90)
    assert np.allclose(vector_crystal_to_cartesian([1,2,3], box),
                       np.array([ 3, 8, 15]))

    box = am.Box(a=3, b=4, c=5, alpha=90, beta=90, gamma=80)
    assert np.allclose(vector_crystal_to_cartesian([1,1,1], box),
                       np.array([[3.69459271, 3.93923101, 5.        ]]))
    # Test 4 indices
    box = am.Box(a=3, b=3, c=5, alpha=90, beta=90, gamma=120)
    assert np.allclose(vector_crystal_to_cartesian([-2, 1, 1, 1], box),
                       np.array([-9, 0, 5]))
    
    box = am.Box(a=3, b=4, c=5, alpha=90, beta=90, gamma=120)
    with pytest.raises(ValueError):
        vector_crystal_to_cartesian([-2, 1, 1, 1], box)

def test_vector_primitive_to_conventional():
    prim = np.array([[ 0.,  1.,  2.],
                     [-3., -2.,  0.],
                     [-3.,  1., -1.],
                     [ 2.,  0., -3.]])
    iset = np.array([[-1.5, -0.5,  0.5],
                     [-0.5, -2.5, -0.5],
                     [-1.5, -0.5, -2.5],
                     [ 2.5,  2.5, -0.5]])
    fset = np.array([[ 0.5,  1. ,  1.5],
                     [-2.5, -1.5, -1. ],
                     [-1. , -2. ,  0. ],
                     [ 1. , -0.5, -1.5]])
    assert np.allclose(vector_primitive_to_conventional(prim, 'i'), iset)
    assert np.allclose(vector_primitive_to_conventional(prim, 'f'), fset)

def test_vector_conventional_to_primitive():
    prim = np.array([[ 0.,  1.,  2.],
                     [-3., -2.,  0.],
                     [-3.,  1., -1.],
                     [ 2.,  0., -3.]])
    iset = np.array([[-1.5, -0.5,  0.5],
                     [-0.5, -2.5, -0.5],
                     [-1.5, -0.5, -2.5],
                     [ 2.5,  2.5, -0.5]])
    fset = np.array([[ 0.5,  1. ,  1.5],
                     [-2.5, -1.5, -1. ],
                     [-1. , -2. ,  0. ],
                     [ 1. , -0.5, -1.5]])
    assert np.allclose(vector_conventional_to_primitive(iset, 'i'), prim)
    assert np.allclose(vector_conventional_to_primitive(fset, 'f'), prim)                     