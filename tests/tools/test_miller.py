# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am
from atomman.tools.miller import *

class TestMiller():

    @property
    def cubbox(self):
        """A cubic box for testing conversions"""
        return am.Box.cubic(a=4)

    @property
    def hexbox(self):
        """A hexagonal box for testing conversions"""
        return am.Box.hexagonal(a=3, c=5)

    @property
    def vector3(self):
        """3 index vectors corresponding to the vector4 4 index vectors."""
        return np.array([[ 1, 2, 3],
                         [ 0, 3, 0],
                         [-3, 0, 1],
                         [-2, 2, 2]])
    
    @property
    def vector4(self):
        """4 index vectors corresponding to the vector3 3 index vectors."""
        return np.array([[ 0, 1,-1, 3],
                         [-1, 2,-1, 0],
                         [-2, 1, 1, 1],
                         [-2, 2, 0, 2]])

    @property
    def plane3(self):
        """3 index planes corresponding to the plane4 4 index planes."""
        return np.array([[ 3, 3, 3],
                         [ 1, 1, 0],
                         [-2, 4, 0],
                         [-1,-3, 5],
                         [ 3, 3,-5]])
    
    @property
    def plane4(self):
        """4 index planes corresponding to the plane3 3 index planes."""
        return np.array([[ 3, 3,-6, 3],
                         [ 1, 1,-2, 0],
                         [-2, 4,-2, 0],
                         [-1,-3, 4, 5],
                         [ 3, 3,-6,-5]])

    @property
    def prim(self):
        """A set of primitive lattice vectors corresponding to iset and fset"""
        return np.array([[ 0.,  1.,  2.],
                         [-3., -2.,  0.],
                         [-3.,  1., -1.],
                         [ 2.,  0., -3.]])

    @property
    def iset(self):
        """A set of body-centered lattice vectors corresponding to prim"""
        return np.array([[-1.5, -0.5,  0.5],
                         [-0.5, -2.5, -0.5],
                         [-1.5, -0.5, -2.5],
                         [ 2.5,  2.5, -0.5]])

    @property
    def fset(self):
        """A set of face-centered lattice vectors corresponding to prim"""
        return np.array([[ 1. ,  0.5,  1.5],
                         [-1.5, -2.5, -1. ],
                         [-2. , -1. ,  0. ],
                         [-0.5,  1. , -1.5]])

    def test_vector3to4(self):
        """test vector3to4 by checking that vector3 -> vector4"""
        # Test giving single vectors
        for vector3, vector4 in zip(self.vector3, self.vector4):
            assert np.all(vector3to4(vector3) == vector4)

        # Test giving array of vectors
        assert np.all(vector3to4(self.vector3) == self.vector4)

        # Test for invalid array shape
        with pytest.raises(ValueError):
            vector3to4([ 1, 2, 3, 4])

    def test_vector4to3(self):
        """test test_vector4to3 by checking that vector4 -> vector3"""
        # Test giving single vectors
        for vector4, vector3 in zip(self.vector4, self.vector3):
            assert np.all(vector4to3(vector4) == vector3)

        # Test giving array of vectors
        assert np.all(vector4to3(self.vector4) == self.vector3)

        # Test for invalid array shape
        with pytest.raises(ValueError):
            vector4to3([ 1, 2, 3])

    def test_plane3to4(self):
        """test plane3to4 by checking that plane3 -> plane4"""
        # Test giving single plane
        for plane3, plane4 in zip(self.plane3, self.plane4):
            assert np.allclose(plane3to4(plane3), plane4)

        # Test giving array of planes
        assert np.allclose(plane3to4(self.plane3), self.plane4)

        # Test for invalid array shape
        with pytest.raises(ValueError):
            plane3to4([ 1, 2, 3, 4])

    def test_plane4to3(self):
        """test plane4to3 by checking that plane4 -> plane3"""
        # Test giving single plane
        for plane4, plane3 in zip(self.plane4, self.plane3):
            assert np.all(plane4to3(plane4) == plane3)

        # Test giving array of planes
        assert np.all(plane4to3(self.plane4) == self.plane3)

        # Test for invalid array shape
        with pytest.raises(ValueError):
            plane4to3([ 1, 2, 3])

    
    
    def test_vector_crystal_to_cartesian(self):
        """
        test vector_crystal_to_cartesian using cubbox and hexbox boxes
        and vector3 and vector4 crystal vector sets.        
        """
        # Define the Cartesian answers
        cubcart = 4 * self.vector3
        hexcart = np.array([[ 0.0, 3*3**0.5,  15.0],
                            [-4.5, 4.5*3**0.5, 0.0],
                            [-9.0, 0.0,        5.0],
                            [-9.0, 3*3**0.5,  10.0]])

        # Test cubic with individual vector3 values
        for vector3, cart in zip(self.vector3, cubcart):
            assert np.allclose(vector_crystal_to_cartesian(vector3, self.cubbox), cart)

        # Test cubic with all vector3 values
        assert np.allclose(vector_crystal_to_cartesian(self.vector3, self.cubbox), cubcart)

        # Test hexagonal with individual vector3 values
        for vector3, cart in zip(self.vector3, hexcart):
            assert np.allclose(vector_crystal_to_cartesian(vector3, self.hexbox), cart)

        # Test hexagonal with all vector3 values   
        assert np.allclose(vector_crystal_to_cartesian(self.vector3, self.hexbox), hexcart)

        # Test hexagonal with vector4
        assert np.allclose(vector_crystal_to_cartesian(self.vector4, self.hexbox), hexcart)
        
        # Check that vector4 fails for cubic
        with pytest.raises(ValueError):
            vector_crystal_to_cartesian(self.vector4, self.cubbox)

    def test_plane_crystal_to_cartesian(self):
        """
        test plane_crystal_to_cartesian using cubbox and hexbox boxes
        and plane3 and plane4 crystal vector sets.  
        """
        # Define the Cartesian answers
        cubcart = (self.plane3.T / np.linalg.norm(self.plane3, axis=1)).T
        hexcart = np.array([[ 0.47891314,  0.8295019 ,  0.28734789],
                            [ 0.5       ,  0.8660254 , -0.        ],
                            [-0.5       ,  0.8660254 , -0.        ],
                            [-0.19487094, -0.78756153,  0.58461282],
                            [ 0.4472136 ,  0.77459667, -0.4472136 ]])

        # Test cubic with individual plane3 values
        for plane3, cart in zip(self.plane3, cubcart):
            assert np.allclose(plane_crystal_to_cartesian(plane3, self.cubbox), cart)

        # Test cubic with all plane3 values
        assert np.allclose(plane_crystal_to_cartesian(self.plane3, self.cubbox), cubcart)

        # Test hexagonal with individual plane3 values
        for plane3, cart in zip(self.plane3, hexcart):
            assert np.allclose(plane_crystal_to_cartesian(plane3, self.hexbox), cart)

        # Test hexagonal with all plane3 values
        assert np.allclose(plane_crystal_to_cartesian(self.plane3, self.hexbox), hexcart)

        # Test hexagonal with plane4
        assert np.allclose(plane_crystal_to_cartesian(self.plane4, self.hexbox), hexcart)
        
        # Check that plane4 fails for cubic
        with pytest.raises(ValueError):
            plane_crystal_to_cartesian(self.plane4, self.cubbox)



    def test_vector_primitive_to_conventional(self):
        """
        test vector_primitive_to_conventional by checking that prim converts
        to prim, iset, or fset based on basis
        """
        assert np.allclose(vector_primitive_to_conventional(self.prim, 'p'), self.prim)
        assert np.allclose(vector_primitive_to_conventional(self.prim, 'i'), self.iset)
        assert np.allclose(vector_primitive_to_conventional(self.prim, 'f'), self.fset)

    def test_vector_conventional_to_primitive(self):
        """
        test vector_conventional_to_primitive by checking that prim, iset,
        and fset convert to prim with the appropriate basis
        """
        assert np.allclose(vector_conventional_to_primitive(self.prim, 'p'), self.prim)
        assert np.allclose(vector_conventional_to_primitive(self.iset, 'i'), self.prim)
        assert np.allclose(vector_conventional_to_primitive(self.fset, 'f'), self.prim)                     