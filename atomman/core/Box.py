# coding: utf-8

# Standard Python libraries
from __future__ import annotations
from copy import deepcopy
import io
from typing import Any, Optional, Union, Tuple

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# local imports
import atomman.unitconvert as uc
from ..tools import vect_angle, miller
from ..region import Shape, Plane

class Box(Shape, object):
    """
    A representation of a triclinic (parallelepiped) box.
    """

    def __init__(self, **kwargs):
        """
        Initializes a Box.  If parameters besides origin are given they
        must completely define the box.  Allowed parameter sets are:

        - no parameters -> box is set to square unit box with origin = [0,0,0].

        - origin. -> Only origin is changed (same as setting origin directly).

        - vects, (and origin).

        - avect, bvect, cvect, (and origin).

        - a, b, c, (alpha, beta, gamma, and origin).

        - lx, ly, lz, (xy, xz, yz, and origin).

        - xlo, xhi, ylo, yhi, zlo, zhi, (xy, xz, and yz).

        - model

        See the description of class methods and attributes for more details
        on the allowed parameters.

        """
        self.__vects = np.eye(3, dtype='float64')
        self.__origin = np.zeros(3, dtype='float64')
        self.__reciprocal_vects = None

        if len(kwargs) > 0:
            if 'model' in kwargs:
                if len(kwargs) > 1:
                    raise ValueError('model cannot be given with other parameters')
                self.model(kwargs['model'])
            else:
                self.set(**kwargs)

    @classmethod
    def cubic(cls, a: float) -> Box:
        """
        Initializes a Box in standard cubic setting using only cubic lattice
        parameters.

        a = b = c, alpha = beta = gamma = 90

        Parameters
        ----------
        - a : float
            The a lattice constant

        Returns
        -------
        atomman.Box
        """
        return cls(a=a, b=a, c=a, alpha=90, beta=90, gamma=90)

    @classmethod
    def hexagonal(cls, a: float, c: float) -> Box:
        """
        Initializes a Box in standard hexagonal setting using only hexagonal lattice
        parameters.

        a = b != c, alpha = beta = 90, gamma = 120

        Parameters
        ----------
        - a : float
            The a lattice constant
        - c : float
            The c lattice constant

        Returns
        -------
        atomman.Box
        """
        if a == c:
            raise ValueError('hexagonal lattice constants must be different')

        return cls(a=a, b=a, c=c, alpha=90, beta=90, gamma=120)

    @classmethod
    def tetragonal(cls, a: float, c: float) -> Box:
        """
        Initializes a Box in standard tetragonal setting using only tetragonal lattice
        parameters.

        a = b != c, alpha = beta = gamma = 90

        Parameters
        ----------
        - a : float
            The a lattice constant
        - c : float
            The c lattice constant

        Returns
        -------
        atomman.Box
        """
        if a == c:
            raise ValueError('tetragonal lattice constants must be different')

        return cls(a=a, b=a, c=c, alpha=90, beta=90, gamma=90)

    @classmethod
    def trigonal(cls, a: float, alpha: float) -> Box:
        """
        Initializes a Box in standard trigonal setting using only trigonal lattice
        parameters.

        a = b = c, alpha = beta = gamma < 120

        Parameters
        ----------
        - a : float
            The a lattice constant
        - alpha : float
            The alpha lattice angle in degrees

        Returns
        -------
        atomman.Box
        """
        if alpha >= 120.0:
            raise ValueError('trigonal alpha angle must be less than 120 degrees')

        return cls(a=a, b=a, c=a, alpha=alpha, beta=alpha, gamma=alpha)

    @classmethod
    def orthorhombic(cls, a: float, b: float, c: float) -> Box:
        """
        Initializes a Box in standard orthorhombic setting using only orthorhombic lattice
        parameters.

        a != b != c, alpha = beta = gamma = 90

        Parameters
        ----------
        - a : float
            The a lattice constant
        - b : float
            The b lattice constant
        - c : float
            The c lattice constant

        Returns
        -------
        atomman.Box
        """
        if a == b or a == c:
            raise ValueError('orthorhombic lattice constants must be different')

        return cls(a=a, b=b, c=c, alpha=90, beta=90, gamma=90)

    @classmethod
    def monoclinic(cls, a: float, b: float, c: float, beta: float) -> Box:
        """
        Initializes a Box in standard monoclinic setting using only monoclinic lattice
        parameters.

        a != b != c, alpha = gamma = 90, beta > 90

        Parameters
        ----------
        - a : float
            The a lattice constant
        - b : float
            The b lattice constant
        - c : float
            The c lattice constant
        - beta : float
            The beta lattice angle in degrees

        Returns
        -------
        atomman.Box
        """
        if a == b or a == c:
            raise ValueError('monoclinic lattice constants must be different')
        if beta <= 90.0:
            raise ValueError('monoclinic angle beta must be greater than 90 degrees')

        return cls(a=a, b=b, c=c, alpha=90, beta=beta, gamma=90)

    @classmethod
    def triclinic(cls, a: float, b: float, c: float,
                  alpha: float, beta: float, gamma: float) -> Box:
        """
        Initializes a Box in standard triclinic setting using only triclinic lattice
        parameters.

        a != b != c, alpha != beta != gamma

        Parameters
        ----------
        - a : float
            The a lattice constant
        - b : float
            The b lattice constant
        - c : float
            The c lattice constant
        - alpha : float
            The alpha lattice angle in degrees
        - beta : float
            The beta lattice angle in degrees
        - gamma : float
            The gamma lattice angle in degrees

        Returns
        -------
        atomman.Box
        """
        if a == b or a == c:
            raise ValueError('monoclinic lattice constants must be different')
        if alpha == beta or alpha == gamma:
            raise ValueError('monoclinic lattice angles must be different')

        return cls(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)

    @property
    def vects(self) -> np.ndarray:
        """numpy.ndarray : Array containing all three box vectors.  Can be set directly."""
        return deepcopy(self.__vects)

    @vects.setter
    def vects(self, value: npt.ArrayLike):
        self.__vects[:] = value

        # Zero out near zero terms
        self.__vects[np.isclose(self.__vects/abs(self.__vects).max(), 0.0, atol=1e-9)] = 0.0

        # Reset reciprocal_vects
        self.__reciprocal_vects = None

    @property
    def reciprocal_vects(self) -> np.ndarray:
        """
        numpy.ndarray : Array of the crystallographic reciprocal box vectors.
        These have not been scaled by the factor of 2 pi.
        """
        if self.__reciprocal_vects is None:
            self.__reciprocal_vects = np.linalg.inv(self.vects).T
        
        return self.__reciprocal_vects

    @property
    def origin(self) -> np.ndarray:
        """numpy.ndarray : Box origin position where vects are added to define the box.  Can be set directly."""
        return deepcopy(self.__origin)

    @origin.setter
    def origin(self, value: npt.ArrayLike):
        self.__origin[:] = value

    @property
    def avect(self) -> np.ndarray:
        """numpy.ndarray : Vector associated with the a box dimension."""
        return self.vects[0]

    @property
    def bvect(self) -> np.ndarray:
        """numpy.ndarray : Vector associated with the b box dimension."""
        return self.vects[1]

    @property
    def cvect(self) -> np.ndarray:
        """numpy.ndarray : Vector associated with the c box dimension."""
        return self.vects[2]

    @property
    def a(self) -> float:
        """float : The a lattice parameter (magnitude of avect)."""
        return (self.__vects[0,0]**2 + self.__vects[0,1]**2 + self.__vects[0,2]**2)**0.5

    @property
    def b(self) -> float:
        """float : The b lattice parameter (magnitude of avect)."""
        return (self.__vects[1,0]**2 + self.__vects[1,1]**2 + self.__vects[1,2]**2)**0.5

    @property
    def c(self) -> float:
        """float : The c lattice parameter (magnitude of avect)."""
        return (self.__vects[2,0]**2 + self.__vects[2,1]**2 + self.__vects[2,2]**2)**0.5

    @property
    def alpha(self) -> float:
        """float : The alpha lattice angle in degrees (angle between bvect and cvect)."""
        return vect_angle(self.__vects[1], self.__vects[2])

    @property
    def beta(self) -> float:
        """float : The beta lattice angle in degrees (angle between avect and cvect)."""
        return vect_angle(self.__vects[0], self.__vects[2])

    @property
    def gamma(self) -> float:
        """float : The gamma lattice angle in degrees (angle between avect and bvect)."""
        return vect_angle(self.__vects[0], self.__vects[1])

    @property
    def lx(self) -> float:
        """float : LAMMPS lx box length (avect[0] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[0,0]

    @property
    def ly(self) -> float:
        """float : LAMMPS ly box length (bvect[1] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[1,1]

    @property
    def lz(self) -> float:
        """float : LAMMPS lz box length (cvect[2] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,2]

    @property
    def xy(self) -> float:
        """float : LAMMPS xy box tilt factor (bvect[0] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[1,0]

    @property
    def xz(self) -> float:
        """float : LAMMPS xz box tilt factor (cvect[0] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,0]

    @property
    def yz(self) -> float:
        """float : LAMMPS yz box tilt factor (cvect[1] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__vects[2,1]

    @property
    def xlo(self) -> float:
        """float : LAMMPS xlo box lo term (origin[0] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[0]

    @property
    def ylo(self) -> float:
        """float : LAMMPS ylo box lo term (origin[1] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[1]

    @property
    def zlo(self) -> float:
        """float : LAMMPS zlo box lo term (origin[2] for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[2]

    @property
    def xhi(self) -> float:
        """float : LAMMPS xhi box hi term (origin[0] + lx for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[0] + self.__vects[0,0]

    @property
    def yhi(self) -> float:
        """float : LAMMPS yhi box hi term (origin[1] + ly for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[1] + self.__vects[1,1]

    @property
    def zhi(self) -> float:
        """float : LAMMPS zhi box hi term (origin[2] + lz for normalized boxes)."""
        assert self.is_lammps_norm(), 'Box is not normalized for LAMMPS style parameters'
        return self.__origin[2] + self.__vects[2,2]

    @property
    def volume(self) -> float:
        """float : The volume of the box."""
        return np.abs(np.dot(self.avect, np.cross(self.bvect, self.cvect)))

    @property
    def planes(self) -> Tuple[Plane]:
        """tuple : The box's planes represented as atomman.region.Plane objects."""
        return (Plane(np.cross(self.cvect, self.bvect), self.origin),
                Plane(np.cross(self.avect, self.cvect), self.origin),
                Plane(np.cross(self.bvect, self.avect), self.origin),
                Plane(np.cross(self.bvect, self.cvect), self.origin + self.avect),
                Plane(np.cross(self.cvect, self.avect), self.origin + self.bvect),
                Plane(np.cross(self.avect, self.bvect), self.origin + self.cvect))

    def __str__(self) -> str:
        """
        The string representation of the box.  Lists the three vectors and origin.
        """
        return '\n'.join(['avect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[0,0], self.__vects[0,1], self.__vects[0,2]),
                          'bvect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[1,0], self.__vects[1,1], self.__vects[1,2]),
                          'cvect =  [%6.3f, %6.3f, %6.3f]' % (self.__vects[2,0], self.__vects[2,1], self.__vects[2,2]),
                          'origin = [%6.3f, %6.3f, %6.3f]' % (self.__origin[0],  self.__origin[1],  self.__origin[2])])

    def model(self,
              model: Union[str, io.IOBase, DM, None] = None,
              length_unit: str = 'angstrom'
              ) -> Optional[DM]:
        """
        Reads or generates a data model for the box.

        Parameters
        ----------
        model : str or DataModelDict, optional
            JSON/XML formatted data, or path to file containing said data.  If
            not given, then a model for the current box will be returned.
        length_unit : str, optional
            Unit of length to save box values in if data model is to be
            generated.  Default value is 'angstrom'.

        Returns
        -------
        DataModelDict.DataModelDict
            A JSON/XML equivalent data model for the box.  Returned if model
            is not given
        """
        # Set values if model given
        if model is not None:

            # Find box element
            model = DM(model).find('box')
            avect = uc.value_unit(model['avect'])
            bvect = uc.value_unit(model['bvect'])
            cvect = uc.value_unit(model['cvect'])
            origin = uc.value_unit(model['origin'])
            self.set(avect=avect, bvect=bvect, cvect=cvect, origin=origin)

        # Return DataModelDict if model not given
        else:
            model = DM()
            model['box'] = DM()
            model['box']['avect'] = uc.model(self.avect, length_unit)
            model['box']['bvect'] = uc.model(self.bvect, length_unit)
            model['box']['cvect'] = uc.model(self.cvect, length_unit)
            model['box']['origin']= uc.model(self.origin, length_unit)

            return model

    def set(self, **kwargs):
        """
        Sets a Box's dimensions.  If parameters besides origin are given they
        must completely define the box.  Allowed parameter sets are:

        - no parameters -> box is set to square unit box with origin = [0,0,0].

        - origin. -> Only origin is changed (same as setting origin directly).

        - vects, (and origin).

        - avect, bvect, cvect, (and origin).

        - a, b, c, (alpha, beta, gamma, and origin).

        - lx, ly, lz, (xy, xz, yz, and origin).

        - xlo, xhi, ylo, yhi, zlo, zhi, (xy, xz, and yz).

        See the description of class methods and attributes for more details
        on the allowed parameters.
        """

        # Set default values if no kwargs given
        if len(kwargs) == 0:
            self.vects = np.eye(3)
            self.origin = np.zeros(3)

        # Set directly if vects given
        elif 'vects' in kwargs:
            vects = kwargs.pop('vects')
            origin = kwargs.pop('origin', [0.0, 0.0, 0.0])
            assert len(kwargs) == 0, 'Invalid arguments'
            self.vects = vects
            self.origin = origin

        # Call set_vectors if vect inputs given
        elif 'avect' in kwargs:
            self.set_vectors(**kwargs)

        # Call set_lengths if length inputs given
        elif 'lx' in kwargs:
            self.set_lengths(**kwargs)

        # Call set_hi_los if hi/lo inputs given
        elif 'xlo' in kwargs:
            self.set_hi_los(**kwargs)

        # Call set_abc if vector magnitudes are given
        elif 'a' in kwargs:
            self.set_abc(**kwargs)

        # Set only origin if given alone
        elif 'origin' in kwargs:
            origin = kwargs.pop('origin')
            assert len(kwargs) == 0, 'Invalid arguments'
            self.origin = origin

        else:
            raise TypeError('Invalid arguments')

    def set_vectors(self,
                    avect: npt.ArrayLike,
                    bvect: npt.ArrayLike,
                    cvect: npt.ArrayLike,
                    origin: Optional[npt.ArrayLike] = None):
        """
        Set the box using the three box vectors.

        Parameters
        ----------
        avect : array-like object
            The 3D vector for the a box dimension.
        bvect : array-like object
            The 3D vector for the b box dimension.
        cvect : array-like object
            The 3D vector for the c box dimension.
        origin : array-like object, optional
            The 3D vector for the box origin position.  Default value is
            (0,0,0).
        """
        # Set default origin
        if origin is None:
            origin = [0.0, 0.0, 0.0]

        # Combine avect, bvect and cvect into vects and set directly
        self.vects = [avect, bvect, cvect]
        self.origin = origin

    def set_abc(self,
                a: float, b: float, c: float,
                alpha: float = 90.0,
                beta: float = 90.0,
                gamma: float = 90.0,
                origin: Optional[npt.ArrayLike] = None):
        """
        Set the box using crystal cell lattice parameters and angles.

        Parameters
        ----------
        a : float
            The a lattice parameter.
        b : float
            The b lattice parameter.
        c : float
            The c lattice parameter.
        alpha : float, optional
            The alpha lattice angle in degrees (angle between b and c vectors).
            Default value is 90.0.
        beta : float, optional
            The beta lattice angle in degrees (angle between a and c vectors).
            Default value is 90.0.
        gamma : float, optional
            The gamma lattice angle in degrees (angle between a and b vectors).
            Default value is 90.0.
        origin : array-like object, optional
            The 3D vector for the box origin position.  Default value is
            (0,0,0).
        """
        # Check that angles are between 0 and 180
        if alpha <= 0 or alpha >=180 or beta <= 0 or beta >= 180 or gamma <=0 or gamma >=180:
            raise ValueError('lattice angles must be between 0 and 180 degrees')

        # Convert to lx, ly, lz, xy, xz, yz
        lx = a
        xy = b * np.cos(gamma * np.pi / 180)
        xz = c * np.cos(beta * np.pi / 180)
        ly = (b**2 - xy**2)**0.5
        yz = (b * c * np.cos(alpha * np.pi / 180) - xy * xz) / ly
        lz = (c**2 - xz**2 - yz**2)**0.5

        # Call set_lengths
        self.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=origin)

    def set_lengths(self,
                    lx: float, ly: float, lz: float,
                    xy: float = 0.0, xz: float = 0.0, yz: float = 0.0,
                    origin: Optional[npt.ArrayLike] = None):
        """
        Set the box using LAMMPS box lengths and tilt factors.

        Parameters
        ----------
        lx : float
            The LAMMPS box length in the x direction.
        ly : float
            The LAMMPS box length in the y direction.
        lz : float
            The LAMMPS box length in the z direction.
        xy : float, optional
            The LAMMPS box tilt factor in the xy direction.  Default value is
            0.0.
        xz : float, optional
            The LAMMPS box tilt factor in the xz direction.  Default value is
            0.0.
        yz : float, optional
            The LAMMPS box tilt factor in the yz direction.  Default value is
            0.0.
        origin : array-like object, optional
            The 3D vector for the box origin position.  Default value is
            (0,0,0).
        """

        assert lx > 0 and ly > 0 and lz > 0, 'box lengths must be positive'

        # Set default origin
        if origin is None:
            origin = [0.0, 0.0, 0.0]

        # Construct vects array
        self.vects = [[lx, 0.0, 0.0],
                      [xy, ly,  0.0],
                      [xz, yz,  lz]]
        self.origin = origin

    def set_hi_los(self,
                   xlo: float, xhi: float,
                   ylo: float, yhi: float,
                   zlo: float, zhi: float,
                   xy: float = 0.0, xz: float = 0.0, yz: float = 0.0):
        """
        Set the box using LAMMPS box hi's, lo's and tilt factors.

        Parameters
        ----------
        xlo : float
            The LAMMPS xlo box lo term.
        xhi : float
            The LAMMPS xhi box hi term.
        ylo : float
            The LAMMPS ylo box lo term.
        yhi : float
            The LAMMPS yhi box hi term.
        zlo : float
            The LAMMPS zlo box lo term.
        zhi : float
            The LAMMPS zhi box hi term.
        xy : float, optional
            The LAMMPS box tilt factor in the xy direction.  Default value is
            0.0.
        xz : float, optional
            The LAMMPS box tilt factor in the xz direction.  Default value is
            0.0.
        yz : float, optional
            The LAMMPS box tilt factor in the yz direction.  Default value is
            0.0.
        """

        # Convert to hi and lo term to lengths and origin
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo
        origin = [xlo, ylo, zlo]

        # Call set_lengths
        self.set_lengths(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=origin)

    def is_lammps_norm(self) -> bool:
        """
        Tests if box is compatible with LAMMPS.
        Note: large box tilt factors not checked.  The LAMMPS command
        'box tilt large' may be needed to run LAMMPS.
        """
        return (self.__vects[0,1] == 0.0
            and self.__vects[0,2] == 0.0
            and self.__vects[1,2] == 0.0
            and self.__vects[0,0] > 0.0
            and self.__vects[1,1] > 0.0
            and self.__vects[2,2] > 0.0)

    def inside(self,
               pos: npt.ArrayLike,
               inclusive: bool = True) -> np.ndarray:
        """
        Indicates if position(s) are inside the shape.
        
        Parameters
        ----------
        pos : array-like object
            Nx3 array of coordinates. 
        inclusive : bool, optional
            Indicates if points on the shape's boundaries are to be included.
            Default value is True.
        
        Returns
        -------
        numpy.NDArray
            N array of bool values: True if inside shape
        """
        # Retrieve the Box's planes
        planes = self.planes

        # Find all points below each plane, i.e. inside the box
        return ( planes[0].below(pos, inclusive=inclusive)
               & planes[1].below(pos, inclusive=inclusive)
               & planes[2].below(pos, inclusive=inclusive)
               & planes[3].below(pos, inclusive=inclusive)
               & planes[4].below(pos, inclusive=inclusive)
               & planes[5].below(pos, inclusive=inclusive))

    def vector_crystal_to_cartesian(self, indices: npt.ArrayLike) -> np.ndarray:
        """
        Converts crystal indices to Cartesian vectors relative
        to the box's lattice vectors. 
        
        Parameters
        ----------
        indices : array-like object
            (..., 3) array of [uvw] Miller crystallographic indices or 
            (..., 4) array of [uvtw] Miller-Bravais crystallographic indices.
    
        Returns
        -------
        np.ndarray of float
            (..., 3) array of Cartesian vectors.
            
        Raises
        ------
        ValueError
            If indices dimensions are not (..., 3) or (..., 4), or if
            hexagonal indices given with non-hexagonal box.
        """
        return miller.vector_crystal_to_cartesian(indices, self)

    def plane_crystal_to_cartesian(self, indices: npt.ArrayLike) -> np.ndarray:
        """
        Converts crystal planar indices to Cartesian plane normal vectors based
        on the box's lattice vectors.  Note: the algorithm used requires that the
        planar indices be integers.
        
        Parameters
        ----------
        indices : array-like object
            (..., 3) array of [hkl] Miller crystallographic indices or 
            (..., 4) array of [hkil] Miller-Bravais crystallographic indices.
        box : atomman.Box
            Box that defines the lattice cell vectors to use. 
    
        Returns
        -------
        np.ndarray of float
            (..., 3) array of Cartesian vectors corresponding to plane normals.
            
        Raises
        ------
        ValueError
            If indices dimensions are not (..., 3) or (..., 4), or if
            hexagonal indices given with non-hexagonal box.
        """
        return miller.plane_crystal_to_cartesian(indices, self)

    def position_relative_to_cartesian(self, relpos: npt.ArrayLike) -> np.ndarray:
        """
        Converts position vectors from relative box coordinates to absolute
        Cartesian coordinates based on the box's vects and origin.

        Parameters
        ----------
        relpos : array-like object
            (..., 3) array of relative position vectors.

        Returns
        -------
        numpy.ndarray
            (..., 3) array of the absolute Cartesian positions corresponding
            to relpos. 

        Raises
        ------
        ValueError
            If relpos dimensions are not (..., 3).
        """
        # Check/convert relpos
        relpos = np.asarray(relpos, dtype=float)
        if relpos.shape[-1] != 3:
            raise ValueError('Invalid position dimensions')

        # Convert and return
        return relpos.dot(self.vects) + self.origin

    def position_cartesian_to_relative(self, cartpos: npt.ArrayLike) -> np.ndarray:
        """
        Converts position vectors from absolute Cartesian coordinates to
        relative box coordinates based on the box's vects and origin.

        Parameters
        ----------
        cartpos : array-like object
            (..., 3) array of Cartesian position vectors.

        Returns
        -------
        numpy.ndarray
            (..., 3) array of the relative positions corresponding
            to cartpos. 

        Raises
        ------
        ValueError
            If cartpos dimensions are not (..., 3).
        """
        # Check/convert cartpos
        value = np.asarray(cartpos, dtype=float)
        if cartpos.shape[-1] != 3:
            raise ValueError('Invalid position dimensions')

        # Convert and return
        return np.inner((value - self.origin), self.reciprocal_vects)

    def iscubic(self, 
                rtol: float = 1e-05,
                atol: float = 1e-08) -> bool:
        """
        Tests if the box is consistent with a standard cubic cell:
        a = b = c
        alpha = beta = gamma = 90
        
        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for testing box parameters. Default value is 1e-5.
        atol : float, optional
            Absolute tolerance for testing box parameters. Default value is 1e-8.
            
        Returns
        -------
        bool
            True if the box is a standard cubic cell, False otherwise.
        """
        return (np.isclose(self.a, self.b, atol=atol, rtol=rtol)
                and np.isclose(self.a, self.c, atol=atol, rtol=rtol)
                and np.isclose(self.alpha, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.beta, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.gamma, 90.0, atol=atol, rtol=rtol))

    def ishexagonal(self,
                    rtol: float = 1e-05,
                    atol: float = 1e-08) -> bool:
        """
        Tests if the box is consistent with a standard hexagonal cell:
        a = b != c
        alpha = beta = 90
        gamma = 120
        
        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for testing box parameters. Default value is 1e-5.
        atol : float, optional
            Absolute tolerance for testing box parameters. Default value is 1e-8.
            
        Returns
        -------
        bool
            True if the box is a standard hexagonal cell, False otherwise.
        """
        return (np.isclose(self.a, self.b, atol=atol, rtol=rtol)
                and np.isclose(self.alpha, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.beta, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.gamma, 120.0, atol=atol, rtol=rtol))

    def istetragonal(self,
                     rtol: float = 1e-05,
                     atol: float = 1e-08) -> bool:
        """
        Tests if the box is consistent with a standard tetragonal cell:
        a = b != c
        alpha = beta = gamma = 90
        
        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for testing box parameters. Default value is 1e-5.
        atol : float, optional
            Absolute tolerance for testing box parameters. Default value is 1e-8.
            
        Returns
        -------
        bool
            True if the box is a standard tetragonal cell, False otherwise.
        """
        return (np.isclose(self.a, self.b, atol=atol, rtol=rtol)
                and not np.isclose(self.a, self.c, atol=atol, rtol=rtol)
                and np.isclose(self.alpha, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.beta, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.gamma, 90.0, atol=atol, rtol=rtol))

    def isrhombohedral(self,
                       rtol: float = 1e-05,
                       atol: float = 1e-08) -> bool:
        """
        Tests if the box is consistent with a standard rhombohedral cell:
        a = b = c
        alpha = beta = gamma != 90
        
        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for testing box parameters. Default value is 1e-5.
        atol : float, optional
            Absolute tolerance for testing box parameters. Default value is 1e-8.
            
        Returns
        -------
        bool
            True if the box is a standard rhombohedral cell, False otherwise.
        """
        return (np.isclose(self.a, self.b, atol=atol, rtol=rtol)
                and np.isclose(self.a, self.c, atol=atol, rtol=rtol)
                and np.isclose(self.alpha, self.beta, atol=atol, rtol=rtol)
                and np.isclose(self.alpha, self.gamma, atol=atol, rtol=rtol)
                and not np.isclose(self.alpha, 90.0, atol=atol, rtol=rtol))

    def isorthorhombic(self,
                       rtol: float = 1e-05,
                       atol: float = 1e-08) -> bool:
        """
        Tests if the box is consistent with a standard orthorhombic cell:
        a != b != c
        alpha = beta = gamma = 90
        
        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for testing box parameters. Default value is 1e-5.
        atol : float, optional
            Absolute tolerance for testing box parameters. Default value is 1e-8.
            
        Returns
        -------
        bool
            True if the box is a standard orthorhombic cell, False otherwise.
        """
        return (not np.isclose(self.a, self.b, atol=atol, rtol=rtol)
                and not np.isclose(self.a, self.c, atol=atol, rtol=rtol)
                and np.isclose(self.alpha, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.beta, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.gamma, 90.0, atol=atol, rtol=rtol))

    def ismonoclinic(self,
                     rtol: float = 1e-05,
                     atol: float = 1e-08) -> bool:
        """
        Tests if the box is consistent with a standard monoclinic cell:
        a != b != c
        alpha = gamma = 90
        beta != 90
        
        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for testing box parameters. Default value is 1e-5.
        atol : float, optional
            Absolute tolerance for testing box parameters. Default value is 1e-8.
            
        Returns
        -------
        bool
            True if box is a standard monoclinic cell, False otherwise.
        """
        return (not np.isclose(self.a, self.b, atol=atol, rtol=rtol)
                and not np.isclose(self.a, self.c, atol=atol, rtol=rtol)
                and np.isclose(self.alpha, 90.0, atol=atol, rtol=rtol)
                and not np.isclose(self.beta, 90.0, atol=atol, rtol=rtol)
                and np.isclose(self.gamma, 90.0, atol=atol, rtol=rtol))

    def istriclinic(self,
                    rtol: float = 1e-05,
                    atol: float = 1e-08) -> bool:
        """
        Tests if the box is consistent with a standard triclinic cell:
        a != b != c
        alpha != 90
        beta != 90
        gamma != 90
        
        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for testing box parameters. Default value is 1e-5.
        atol : float, optional
            Absolute tolerance for testing box parameters. Default value is 1e-8.
            
        Returns
        -------
        bool
            True if box is a standard triclinic cell, False otherwise.
        """
        return (not np.isclose(self.a, self.b, atol=atol, rtol=rtol)
                and not np.isclose(self.a, self.c, atol=atol, rtol=rtol)
                and not np.isclose(self.alpha, self.beta, atol=atol, rtol=rtol)
                and not np.isclose(self.alpha, self.gamma, atol=atol, rtol=rtol))

    def identifyfamily(self, 
                       rtol: float = 1e-05,
                       atol: float = 1e-08) -> Optional[str]:
        """
        Tests if the box is consistent with a standard representation
        of a crystal system cell.
        
        Parameters
        ----------
        rtol : float, optional
            Relative tolerance for testing box parameters. Default value is 1e-5.
        atol : float, optional
            Absolute tolerance for testing box parameters. Default value is 1e-8.
            
        Returns
        -------
        str or None
            'cubic', 'hexagonal', 'tetragonal', 'rhombohedral', 'orthorhombic',
            'monoclinic' or 'triclinic' if it matches any. None if no matches.
            
        Raises
        ------
        ValueError
            If the box is not consistent with a standard cell.
        """
        if self.iscubic(rtol=rtol, atol=atol):
            return 'cubic'
        elif self.ishexagonal(rtol=rtol, atol=atol):
            return 'hexagonal'
        elif self.istetragonal(rtol=rtol, atol=atol):
            return 'tetragonal'
        elif self.isrhombohedral(rtol=rtol, atol=atol):
            return 'rhombohedral'
        elif self.isorthorhombic(rtol=rtol, atol=atol):
            return 'orthorhombic'
        elif self.ismonoclinic(rtol=rtol, atol=atol):
            return 'monoclinic'
        elif self.istriclinic(rtol=rtol, atol=atol):
            return 'triclinic'
        else:
            None

    def d_hkl(self, plane_hkl):
        """
        Compute the interplanar spacing for a lattice plane.  Note this feature
        is expected to be used for only unit cell boxes!

        Parameters
        ----------
        plane_hkl : array-like
            The 3 index Miller or 4 index Miller-Bravais lattice plane.  Values
            typically are ints, but can be floats if the cell basis is not
            primitive.

        Returns
        -------
        d_hkl : float
            The spacing between the lattice planes of the given type.
        """
        # Check hkl values
        plane_hkl = np.asarray(plane_hkl)
        if plane_hkl.shape == (4,):
            if not self.ishexagonal():
                raise ValueError('Box is not hexagonal: cannot use 4 index Miller-Bravais plane')
            plane_hkl = miller.plane4to3(plane_hkl)
        elif plane_hkl.shape != (3,):
            raise ValueError('plane_hkl must have 3 or 4 indices')
        h = plane_hkl[0]
        k = plane_hkl[1]
        l = plane_hkl[2]

        a = self.a
        b = self.b
        c = self.c

        # Compute sin and cos of the angles
        sin_alpha = np.sin(self.alpha / 180 * np.pi)
        sin_beta =  np.sin(self.beta  / 180 * np.pi)
        sin_gamma = np.sin(self.gamma / 180 * np.pi)
        cos_alpha = np.cos(self.alpha / 180 * np.pi)
        cos_beta =  np.cos(self.beta  / 180 * np.pi)
        cos_gamma = np.cos(self.gamma / 180 * np.pi)

        # Compute the triclinic 1/d_hkl^2 formula
        numerator = (
            h**2 / a**2 * sin_alpha**2 +
            k**2 / b**2 * sin_beta**2 +
            l**2 / c**2 * sin_gamma**2 +
            (2 * k * l / b * c) * (cos_beta * cos_gamma - cos_alpha) +
            (2 * h * l / a * c) * (cos_gamma * cos_alpha - cos_beta) +
            (2 * h * k / a * b) * (cos_alpha * cos_beta - cos_gamma))
        denominator = (1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2 +
                       2 * cos_alpha * cos_beta * cos_gamma)
        
        return (numerator / denominator)**-0.5