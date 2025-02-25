# coding: utf-8
# Standard Python libraries
from __future__ import annotations
from copy import deepcopy
import io
from typing import Optional, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
from ..tools import axes_check
import atomman.unitconvert as uc

class ElasticConstants(object):
    """Class for storing and converting elastic constant values"""

    def __init__(self, **kwargs):
        """
        Initializes an ElasticConstants instance from one of the parameter
        options.
        
        Parameters
        ----------
        Cij : numpy.ndarray
            (6, 6) array of Voigt representation of elastic stiffness.
        Sij : numpy.ndarray
            (6, 6) array of Voigt representation of elastic compliance.
        Cij9 : numpy.ndarray
            (9, 9) array representation of elastic stiffness.
        Cijkl : numpy.ndarray
            (3, 3, 3, 3) array representation of elastic stiffness.
        Sijkl : numpy.ndarray
            (3, 3, 3, 3) array representation of elastic compliance.
        model : DataModelDict, string, or file-like object 
            Data model containing elastic constants.
        C11, C12, ... C66 : float
            Individual components of Cij for a standardized representation:
            isotropic: C11, C12, C44 (2*C44=C11-C12)
            cubic: C11, C12, C44
            hexagonal: C11, C12, C13, C33, C44, C66 (2*C66=C11-C12)
            tetragonal: C11, C12, C13, C16, C33, C44, C66 (C16 optional)
            rhombohedral: C11, C12, C13, C14, C15, C33, C44, C66 (2*C66=C11-C12, C15 optional)
            orthorhombic: C11, C12, C13, C22, C23, C33, C44, C55, C66
            monoclinic: C11, C12, C13, C15, C22, C23, C25, C33, C35, C44, C46, C55, C66
            triclinic: all Cij where i <= j
        M, Lame, mu, E, nu, K : float
            Individual isotropic elastic moduli that can be used instead.
        verify : bool, optional
            Used with hexagonal and rhombohedral Cij sets.  If True (default) then values
            of the non-independent moduli C11, C12, C66 will be checked for compatibility
            if all three are given.
        """
        # Pop verify if needed
        verify = kwargs.pop('verify', True)

        # Initialize for no arguments
        if len(kwargs) == 0:
            self.__c_ij = np.zeros((6,6), dtype='float64')

        # Initialize for matrix arguments
        elif 'Cij' in kwargs:
            assert len(kwargs) == 1, 'Cij cannot be specified with other keyword arguments'
            self.Cij = kwargs['Cij']
        elif 'Sij' in kwargs:
            assert len(kwargs) == 1, 'Sij cannot be specified with other keyword arguments'
            self.Sij = kwargs['Sij']
        elif 'Cij9' in kwargs:
            assert len(kwargs) == 1, 'Cij9 cannot be specified with other keyword arguments'
            self.Cij9 = kwargs['Cij9']
        elif 'Cijkl' in kwargs:
            assert len(kwargs) == 1, 'Cijkl cannot be specified with other keyword arguments'
            self.Cijkl = kwargs['Cijkl']
        elif 'Sijkl' in kwargs:
            assert len(kwargs) == 1, 'Sijkl cannot be specified with other keyword arguments'
            self.Sijkl = kwargs['Sijkl']

        # Initialize using data model
        elif 'model' in kwargs:
            self.model(**kwargs)

        # Initialize for individually specified parameters by standard representation
        elif ('C24' in kwargs or 'C26' in kwargs or 'C34' in kwargs or
              'C36' in kwargs or 'C45' in kwargs or 'C56' in kwargs):
            self.triclinic(**kwargs)
        elif 'C25' in kwargs or 'C35' in kwargs or 'C46' in kwargs:
            self.monoclinic(**kwargs)
        elif 'C22' in kwargs or 'C23' in kwargs or 'C55' in kwargs:
            self.orthorhombic(**kwargs)
        elif 'C14' in kwargs:
            self.rhombohedral(verify=verify, **kwargs)
        elif 'C13' in kwargs or 'C33' in kwargs or 'C66' in kwargs:
            if 'C16' in kwargs or ('C11' in kwargs and 'C12' in kwargs and 'C66' in kwargs):
                self.tetragonal(**kwargs)
            else:
                self.hexagonal(verify=verify, **kwargs)
        elif len(kwargs) == 3:
            self.cubic(**kwargs)
        elif len(kwargs) == 2:
            self.isotropic(**kwargs)
        else:
            raise TypeError('Invalid argument keywords')

    def __str__(self)  -> str:
        """Calling string returns str(self.Cij)."""
        return str(self.Cij)

    @property
    def Cij(self) -> np.ndarray:
        """The stiffness constants in Voigt 6x6 format"""
        return deepcopy(self.__c_ij)

    @Cij.setter
    def Cij(self, value: npt.ArrayLike):
        value = np.asarray(value, dtype='float64')
        assert value.shape == (6,6),  'Cij must be 6x6'

        # Zero out near-zero terms
        assert value.max() > 0.0, 'Cij values not valid'
        value[np.isclose(value/value.max(), 0.0, atol=1e-9)] = 0.0

        # Check symmetry
        for i in range(6):
            for j in range(i):
                assert np.isclose(value[i,j], value[j,i], atol=1e-9), '6x6 matrix not symmetric!'
        self.__c_ij = value

    @property
    def Sij(self) -> np.ndarray:
        """The compliance constants in Voigt 6x6 format"""
        return np.linalg.inv(self.Cij)

    @Sij.setter
    def Sij(self, value: npt.ArrayLike):
        value = np.asarray(value, dtype='float64')
        assert value.shape == (6,6),  'Sij must be 6x6'
        self.Cij = np.linalg.inv(value)

    @property
    def Cij9(self) -> np.ndarray:
        """The stiffness constants in 9x9 format"""
        c = self.Cij
        return np.array([[c[0,0],c[0,1],c[0,2],c[0,3],c[0,4],c[0,5],c[0,3],c[0,4],c[0,5]],
                         [c[1,0],c[1,1],c[1,2],c[1,3],c[1,4],c[1,5],c[1,3],c[1,4],c[1,5]],
                         [c[2,0],c[2,1],c[2,2],c[2,3],c[2,4],c[2,5],c[2,3],c[2,4],c[2,5]],
                         [c[3,0],c[3,1],c[3,2],c[3,3],c[3,4],c[3,5],c[3,3],c[3,4],c[3,5]],
                         [c[4,0],c[4,1],c[4,2],c[4,3],c[4,4],c[4,5],c[4,3],c[4,4],c[4,5]],
                         [c[5,0],c[5,1],c[5,2],c[5,3],c[5,4],c[5,5],c[5,3],c[5,4],c[5,5]],
                         [c[3,0],c[3,1],c[3,2],c[3,3],c[3,4],c[3,5],c[3,3],c[3,4],c[3,5]],
                         [c[4,0],c[4,1],c[4,2],c[4,3],c[4,4],c[4,5],c[4,3],c[4,4],c[4,5]],
                         [c[5,0],c[5,1],c[5,2],c[5,3],c[5,4],c[5,5],c[5,3],c[5,4],c[5,5]]])

    @Cij9.setter
    def Cij9(self, value: npt.ArrayLike):
        value = np.asarray(value, dtype='float64')
        assert value.shape == (9,9), 'Cij9 must be 9x9'

        # Check symmetry
        for i in range(6, 9):
            for j in range(9):
                assert value[i,j] == value[i-3, j]
                assert value[j,i] == value[j, i-3]
        self.Cij = value[:6, :6]

    @property
    def Cijkl(self) -> np.ndarray:
        """The stiffness constants in 3x3x3x3 format"""
        c = self.Cij
        return np.array([[[[c[0,0],c[0,5],c[0,4]], [c[0,5],c[0,1],c[0,3]], [c[0,4],c[0,3],c[0,2]]],
                          [[c[5,0],c[5,5],c[5,4]], [c[5,5],c[5,1],c[5,3]], [c[5,4],c[5,3],c[5,2]]],
                          [[c[4,0],c[4,5],c[4,4]], [c[4,5],c[4,1],c[4,3]], [c[4,4],c[4,3],c[4,2]]]],

                         [[[c[5,0],c[5,5],c[5,4]], [c[5,5],c[5,1],c[5,3]], [c[5,4],c[5,3],c[5,2]]],
                          [[c[1,0],c[1,5],c[1,4]], [c[1,5],c[1,1],c[1,3]], [c[1,4],c[1,3],c[1,2]]],
                          [[c[3,0],c[3,5],c[3,4]], [c[3,5],c[3,1],c[3,3]], [c[3,4],c[3,3],c[3,2]]]],

                         [[[c[4,0],c[4,5],c[4,4]], [c[4,5],c[4,1],c[4,3]], [c[4,4],c[4,3],c[4,2]]],
                          [[c[3,0],c[3,5],c[3,4]], [c[3,5],c[3,1],c[3,3]], [c[3,4],c[3,3],c[3,2]]],
                          [[c[2,0],c[2,5],c[2,4]], [c[2,5],c[2,1],c[2,3]], [c[2,4],c[2,3],c[2,2]]]]])

    @Cijkl.setter
    def Cijkl(self, value: npt.ArrayLike):
        c = np.asarray(value, dtype='float64')
        assert c.shape == (3,3,3,3),  'Cijkl must be 3x3x3x3'
        assert c.max() > 0.0, 'Cij values not valid'
        # Check symmetry
        indexes = np.array([[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]], dtype=int)
        for ij in range(6):
            for kl in range(ij, 6):
                i, j, k, l = indexes[ij,0], indexes[ij,1], indexes[kl,0], indexes[kl,1]
                assert np.isclose(c[i,j,k,l], c[j,i,k,l])
                assert np.isclose(c[i,j,k,l], c[j,i,l,k])
                assert np.isclose(c[i,j,k,l], c[k,l,j,i])
                assert np.isclose(c[i,j,k,l], c[l,k,j,i])
                assert np.isclose(c[i,j,k,l], c[i,j,l,k])
                assert np.isclose(c[i,j,k,l], c[k,l,i,j])
                assert np.isclose(c[i,j,k,l], c[l,k,i,j])

        self.Cij = np.array([[c[0,0,0,0], c[0,0,1,1], c[0,0,2,2], c[0,0,1,2], c[0,0,0,2], c[0,0,0,1]],
                             [c[1,1,0,0], c[1,1,1,1], c[1,1,2,2], c[1,1,1,2], c[1,1,0,2], c[1,1,0,1]],
                             [c[2,2,0,0], c[2,2,1,1], c[2,2,2,2], c[2,2,1,2], c[2,2,0,2], c[2,2,0,1]],
                             [c[1,2,0,0], c[1,2,1,1], c[1,2,2,2], c[1,2,1,2], c[1,2,0,2], c[1,2,0,1]],
                             [c[0,2,0,0], c[0,2,1,1], c[0,2,2,2], c[0,2,1,2], c[0,2,0,2], c[0,2,0,1]],
                             [c[0,1,0,0], c[0,1,1,1], c[0,1,2,2], c[0,1,1,2], c[0,1,0,2], c[0,1,0,1]]])

    @property
    def Sijkl(self) -> np.ndarray:
        """The compliance constants in 3x3x3x3 format"""
        s = self.Sij
        s[3:,:] = s[3:,:]/2.
        s[:,3:] = s[:,3:]/2.
        return np.array([[[[s[0,0],s[0,5],s[0,4]], [s[0,5],s[0,1],s[0,3]], [s[0,4],s[0,3],s[0,2]]],
                          [[s[5,0],s[5,5],s[5,4]], [s[5,5],s[5,1],s[5,3]], [s[5,4],s[5,3],s[5,2]]],
                          [[s[4,0],s[4,5],s[4,4]], [s[4,5],s[4,1],s[4,3]], [s[4,4],s[4,3],s[4,2]]]],

                         [[[s[5,0],s[5,5],s[5,4]], [s[5,5],s[5,1],s[5,3]], [s[5,4],s[5,3],s[5,2]]],
                          [[s[1,0],s[1,5],s[1,4]], [s[1,5],s[1,1],s[1,3]], [s[1,4],s[1,3],s[1,2]]],
                          [[s[3,0],s[3,5],s[3,4]], [s[3,5],s[3,1],s[3,3]], [s[3,4],s[3,3],s[3,2]]]],

                         [[[s[4,0],s[4,5],s[4,4]], [s[4,5],s[4,1],s[4,3]], [s[4,4],s[4,3],s[4,2]]],
                          [[s[3,0],s[3,5],s[3,4]], [s[3,5],s[3,1],s[3,3]], [s[3,4],s[3,3],s[3,2]]],
                          [[s[2,0],s[2,5],s[2,4]], [s[2,5],s[2,1],s[2,3]], [s[2,4],s[2,3],s[2,2]]]]])

    @Sijkl.setter
    def Sijkl(self, value: npt.ArrayLike):
        s = np.asarray(value, dtype='float64')
        assert s.shape == (3,3,3,3),  'Sijkl must be 3x3x3x3'

        # Check symmetry
        indexes = np.array([[0,0], [1,1], [2,2], [1,2], [0,2], [0,1]], dtype=int)
        for ij in range(6):
            for kl in range(ij, 6):
                i, j, k, l = indexes[ij,0], indexes[ij,1], indexes[kl,0], indexes[kl,1]
                assert np.isclose(s[i,j,k,l], s[j,i,k,l])
                assert np.isclose(s[i,j,k,l], s[j,i,l,k])
                assert np.isclose(s[i,j,k,l], s[k,l,j,i])
                assert np.isclose(s[i,j,k,l], s[l,k,j,i])
                assert np.isclose(s[i,j,k,l], s[i,j,l,k])
                assert np.isclose(s[i,j,k,l], s[k,l,i,j])
                assert np.isclose(s[i,j,k,l], s[l,k,i,j])

        self.Sij = np.array([[   s[0,0,0,0],    s[0,0,1,1],    s[0,0,2,2], 2.*s[0,0,1,2], 2.*s[0,0,0,2], 2.*s[0,0,0,1]],
                             [   s[1,1,0,0],    s[1,1,1,1],    s[1,1,2,2], 2.*s[1,1,1,2], 2.*s[1,1,0,2], 2.*s[1,1,0,1]],
                             [   s[2,2,0,0],    s[2,2,1,1],    s[2,2,2,2], 2.*s[2,2,1,2], 2.*s[2,2,0,2], 2.*s[2,2,0,1]],
                             [2.*s[1,2,0,0], 2.*s[1,2,1,1], 2.*s[1,2,2,2], 4.*s[1,2,1,2], 4.*s[1,2,0,2], 4.*s[1,2,0,1]],
                             [2.*s[0,2,0,0], 2.*s[0,2,1,1], 2.*s[0,2,2,2], 4.*s[0,2,1,2], 4.*s[0,2,0,2], 4.*s[0,2,0,1]],
                             [2.*s[0,1,0,0], 2.*s[0,1,1,1], 2.*s[0,1,2,2], 4.*s[0,1,1,2], 4.*s[0,1,0,2], 4.*s[0,1,0,1]]])

    def transform(self,
                  axes: npt.ArrayLike,
                  tol: float = 1e-8) -> ElasticConstants:
        """
        Transforms the elastic constant matrix based on the supplied axes.
        
        Parameters
        ----------
        axes : numpy.ndarray
            (3, 3) array giving three right-handed orthogonal vectors to use
            for transforming.
        tol : float, optional
            Relative tolerance to use in identifying near-zero terms.
            
        Returns
        -------
        ElasticConstants
            A new ElasticConstants object that has been transformed.
        """
        axes = np.asarray(axes, dtype='float64')
        T = axes_check(axes)

        Q = np.einsum('km,ln->mnkl', T, T)
        C = np.einsum('ghij,ghmn,mnkl->ijkl', Q, self.Cijkl, Q)
        C[abs(C / C.max()) < tol] = 0.0

        return ElasticConstants(Cijkl=C)

    def isotropic(self, *,
                  C11: Optional[float] = None,
                  C12: Optional[float] = None,
                  C44: Optional[float] = None,
                  M: Optional[float] = None,
                  Lame: Optional[float] = None,
                  mu: Optional[float] = None,
                  E: Optional[float] = None,
                  nu: Optional[float] = None,
                  K: Optional[float] = None):
        """
        Set values with two independent isotropic moduli.
        
        Parameters
        ----------
        C11 : float, optional
            C11 component of Cij.
        C12 : float, optional
            C12 component of Cij.
        C44 : float, optional
            C44 component of Cij.
        M : float, optional
            P-wave modulus(Equivalent to C11).
        Lame : float, optional
            Lame's first parameter (Equivalent to C12).
        mu : float, optional
            Shear modulus (Equivalent to C44).
        E : float, optional
            Young's modulus
        nu : float, optional
            Poisson's ratio
        K : float, optional
            Bulk modulus
        """
        # Count moduli with values
        kwargcount = sum([C11 is not None, C12 is not None, C44 is not None,
                          M is not None, Lame is not None, mu is not None,
                          E is not None, nu is not None, K is not None])
        if kwargcount != 2:
            raise TypeError('isotropic requires exactly two independent moduli')

        # Handle equivalent terms
        if M is not None:
            if C11 is not None:
                raise ValueError('C11 and M are not independent moduli')
            C11 = M
        if Lame is not None:
            if C12 is not None:
                raise ValueError('C12 and Lame are not independent moduli')
            C12 = Lame
        if mu is not None:
            if C44 is not None:
                raise ValueError('C44 and mu are not independent moduli')
            C44 = mu

        # C11 + (something) combinations
        if C11 is not None:

            # Get C44 from C11 and C12
            if C12 is not None:
                C44 = (C11 - C12) / 2

            # Find C44 then get C12
            else:
                if E is not None:
                    S = (E**2 + 9 * C11**2 - 10 * E * C11)**0.5
                    C44 = (3 * C11 + E - S) / 8
                elif nu is not None:
                    C44 = C11 * (1 - 2 * nu) / (2 * (1 - nu))
                elif K is not None:
                    C44 = 3 * (C11 - K) / 4
                # else: C44 is not None

                C12 = C11 - 2 * C44

        # Combinations without C11: find C12 and C44 then get C11
        else:
            # C12 + (something) combinations
            if C12 is not None:
                if E is not None:
                    R = (E**2 + 9 * C12**2 + 2 * E * C12)**0.5
                    C44 = (E - 3 * C12 + R) / 4
                elif nu is not None:
                    C44 = C12 * (1 - 2 * nu) / (2 * nu)
                elif K is not None:
                    C44 = 3 * (K - C12) / 2
                # else: C44 is not None

            # C44 + (something) combinations
            elif C44 is not None:
                if E is not None:
                    C12 = C44 * (E - 2 * C44) / (3 * C44 - E)
                elif nu is not None:
                    C12 = 2 * C44 * nu / (1 - 2 * nu)
                else: # K is not None
                    C12 = K - 2 * C44 / 3

            # Others + (something) combinations
            elif E is not None:
                if nu is not None:
                    C12 = E * nu / ((1 + nu) * (1 - 2 * nu))
                    C44 = E / (2 * (1 + nu))
                else: # K is not None
                    C12 = 3 * K * (3 * K - E) / (9 * K - E)
                    C44 = 3 * K * E / (9 * K - E)
            elif nu is not None:
                #if K is not None:
                C12 = 3 * K * nu / (1 + nu)
                C44 = 3 * K * (1 - 2 * nu) / (2 * (1 + nu))

            C11 = C12 + 2 * C44

        # Build Cij array
        self.Cij = np.array([[C11, C12, C12, 0.0, 0.0, 0.0],
                             [C12, C11, C12, 0.0, 0.0, 0.0],
                             [C12, C12, C11, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, C44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, C44, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, C44]])

    def cubic(self, *,
              C11: Optional[float] = None,
              C12: Optional[float] = None,
              C44: Optional[float] = None):
        """
        Set values with three independent cubic moduli.
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C44 : float
            C44 component of Cij.
        """
        if C11 is None or C12 is None or C44 is None:
            raise TypeError('cubic style requires C11, C12, and C66')

        # Build Cij array
        self.Cij = np.array([[C11, C12, C12, 0.0, 0.0, 0.0],
                             [C12, C11, C12, 0.0, 0.0, 0.0],
                             [C12, C12, C11, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, C44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, C44, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, C44]])

    def hexagonal(self, *,
                  C11: Optional[float] = None,
                  C12: Optional[float] = None,
                  C13: Optional[float] = None,
                  C33: Optional[float] = None,
                  C44: Optional[float] = None,
                  C66: Optional[float] = None,
                  verify: bool = True):
        """
        Set values with five independent hexagonal moduli.
        (2 * C66 = C11 - C12)
        
        Parameters
        ----------
        C11 : float, optional
            C11 component of Cij.
        C12 : float, optional
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C33 : float
            C33 component of Cij.
        C44 : float
            C44 component of Cij.
        C66 : float, optional
            C66 component of Cij.
        verify : bool, optional
            If True (default), the values of the non-independent moduli C11,
            C12, and C66 will be checked for compatibility if all are given.
        """
        if C13 is None or C33 is None or C44 is None:
            raise TypeError('hexagonal style requires C13, C33, and C44')

        # Calculate missing C11, C12 or C66
        if C11 is None:
            if C12 is None or C66 is None:
                raise TypeError('hexagonal style requires at least two of C11, C12, and C66')
            C11 = 2 * C66 + C12
        elif C12 is None:
            if C66 is None:
                raise TypeError('hexagonal style requires at least two of C11, C12, and C66')
            C12 = C11 - 2 * C66
        elif C66 is None:
            C66 = (C11 - C12) / 2

        # Verify 2 * C66 = C11 - C12
        elif verify:
            if not np.isclose(2 * C66, C11 - C12):
                raise ValueError('dependent values not compatible: C11-C12 != 2*C66')

        # Build Cij array
        self.Cij = np.array([[C11, C12, C13, 0.0, 0.0, 0.0],
                             [C12, C11, C13, 0.0, 0.0, 0.0],
                             [C13, C13, C33, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, C44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, C44, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, C66]])

    def rhombohedral(self, *,
                     C11: Optional[float] = None,
                     C12: Optional[float] = None,
                     C13: Optional[float] = None,
                     C14: Optional[float] = None,
                     C15: float = 0.0,
                     C33: Optional[float] = None,
                     C44: Optional[float] = None,
                     C66: Optional[float] = None,
                     verify: bool = True):
        """
        Set values with six or seven independent rhombohedral moduli. 
        (2 * C66 = C11 - C12)
        
        Parameters
        ----------
        C11 : float, optional
            C11 component of Cij.
        C12 : float, optional
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C14 : float
            C14 component of Cij.
        C15 : float, optional
            C15 component of Cij.
        C33 : float
            C33 component of Cij.
        C44 : float
            C44 component of Cij.
        C66 : float, optional
            C66 component of Cij.
        verify : bool, optional
            If True, the values of the non-independent moduli C11, C12, C66
            will be checked for compatibility.
        """
        if C13 is None or C14 is None or C33 is None or C44 is None:
            raise TypeError('rhombohedral style requires C13, C14, C33, and C44')

        # Calculate missing C11, C12 or C66
        if C11 is None:
            if C12 is None or C66 is None:
                raise TypeError('rhombohedral style requires at least two of C11, C12, and C66')
            C11 = 2 * C66 + C12
        elif C12 is None:
            if C66 is None:
                raise TypeError('rhombohedral style requires at least two of C11, C12, and C66')
            C12 = C11 - 2 * C66
        elif C66 is None:
            C66 = (C11 - C12) / 2

        # Verify 2 * C66 = C11 - C12
        elif verify:
            if not np.isclose(2 * C66, C11 - C12):
                raise ValueError('dependent values not compatible: C11-C12 != 2*C66')

        # Build Cij array
        self.Cij = np.array([[C11, C12, C13, C14, C15, 0.0],
                             [C12, C11, C13,-C14,-C15, 0.0],
                             [C13, C13, C33, 0.0, 0.0, 0.0],
                             [C14,-C14, 0.0, C44, 0.0,-C15],
                             [C15,-C15, 0.0, 0.0, C44, C14],
                             [0.0, 0.0, 0.0,-C15, C14, C66]])

    def tetragonal(self, *,
                   C11: Optional[float] = None,
                   C12: Optional[float] = None,
                   C13: Optional[float] = None,
                   C16: float = 0.0,
                   C33: Optional[float] = None,
                   C44: Optional[float] = None,
                   C66: Optional[float] = None):
        """
        Set values with six or seven independent tetragonal moduli.
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C16 : float, optional
            C16 component of Cij.
        C33 : float
            C33 component of Cij.
        C44 : float
            C44 component of Cij.
        C66 : float
            C66 component of Cij.
        """
        if (C11 is None or C12 is None or C13 is None or
            C33 is None or C44 is None or C66 is None):
            raise TypeError('tetragonal style requires C11, C12, C13, C33, C44, and C66')

        # Build Cij array
        self.Cij = np.array([[C11, C12, C13, 0.0, 0.0, C16],
                             [C12, C11, C13, 0.0, 0.0,-C16],
                             [C13, C13, C33, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, C44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, C44, 0.0],
                             [C16,-C16, 0.0, 0.0, 0.0, C66]])

    def orthorhombic(self, *,
                     C11: Optional[float] = None,
                     C12: Optional[float] = None,
                     C13: Optional[float] = None,
                     C22: Optional[float] = None,
                     C23: Optional[float] = None,
                     C33: Optional[float] = None,
                     C44: Optional[float] = None,
                     C55: Optional[float] = None,
                     C66: Optional[float] = None):
        """
        Set values with nine independent orthorhombic moduli.
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C22 : float
            C22 component of Cij.
        C23 : float
            C23 component of Cij.
        C33 : float
            C33 component of Cij.
        C44 : float
            C44 component of Cij.
        C55 : float
            C55 component of Cij.
        C66 : float
            C66 component of Cij.
        """
        if (C11 is None or C12 is None or C13 is None or
            C22 is None or C23 is None or C33 is None or
            C44 is None or C55 is None or C66 is None):
            raise TypeError('orthorhombic style requires C11, C12, C13, C22, C23, C33, C44, C55, and C66')

        # Build Cij array
        self.Cij = np.array([[C11, C12, C13, 0.0, 0.0, 0.0],
                             [C12, C22, C23, 0.0, 0.0, 0.0],
                             [C13, C23, C33, 0.0, 0.0, 0.0],
                             [0.0, 0.0, 0.0, C44, 0.0, 0.0],
                             [0.0, 0.0, 0.0, 0.0, C55, 0.0],
                             [0.0, 0.0, 0.0, 0.0, 0.0, C66]])

    def monoclinic(self, *,
                   C11: Optional[float] = None,
                   C12: Optional[float] = None,
                   C13: Optional[float] = None,
                   C15: Optional[float] = None,
                   C22: Optional[float] = None,
                   C23: Optional[float] = None,
                   C25: Optional[float] = None,
                   C33: Optional[float] = None,
                   C35: Optional[float] = None,
                   C44: Optional[float] = None,
                   C46: Optional[float] = None,
                   C55: Optional[float] = None,
                   C66: Optional[float] = None):
        """
        Set values with thirteen independent monoclinic moduli.
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C15 : float
            C15 component of Cij.
        C22 : float
            C22 component of Cij.
        C23 : float
            C23 component of Cij.
        C25 : float
            C25 component of Cij.
        C33 : float
            C33 component of Cij.
        C35 : float
            C35 component of Cij.
        C44 : float
            C44 component of Cij.
        C46 : float
            C46 component of Cij.
        C55 : float
            C55 component of Cij.
        C66 : float
            C66 component of Cij.
        """
        if (C11 is None or C12 is None or C13 is None or C15 is None or
            C22 is None or C23 is None or C25 is None or
            C33 is None or C35 is None or C44 is None or
            C46 is None or C55 is None or C66 is None):
            raise TypeError('monoclinic style requires C11, C12, C13, C15, C22, C23, C25, C33, C35, C44, C46, C55, and C66')

        # Build Cij array
        self.Cij = np.array([[C11, C12, C13, 0.0, C15, 0.0],
                             [C12, C22, C23, 0.0, C25, 0.0],
                             [C13, C23, C33, 0.0, C35, 0.0],
                             [0.0, 0.0, 0.0, C44, 0.0, C46],
                             [C15, C25, C35, 0.0, C55, 0.0],
                             [0.0, 0.0, 0.0, C46, 0.0, C66]])

    def triclinic(self, *,
                  C11: Optional[float] = None,
                  C12: Optional[float] = None,
                  C13: Optional[float] = None,
                  C14: Optional[float] = None,
                  C15: Optional[float] = None,
                  C16: Optional[float] = None,
                  C22: Optional[float] = None,
                  C23: Optional[float] = None,
                  C24: Optional[float] = None,
                  C25: Optional[float] = None,
                  C26: Optional[float] = None,
                  C33: Optional[float] = None,
                  C34: Optional[float] = None,
                  C35: Optional[float] = None,
                  C36: Optional[float] = None,
                  C44: Optional[float] = None,
                  C45: Optional[float] = None,
                  C46: Optional[float] = None,
                  C55: Optional[float] = None,
                  C56: Optional[float] = None,
                  C66: Optional[float] = None):
        """
        Set values with twenty one independent triclinic moduli
        
        Parameters
        ----------
        C11 : float
            C11 component of Cij.
        C12 : float
            C12 component of Cij.
        C13 : float
            C13 component of Cij.
        C14 : float
            C14 component of Cij.
        C15 : float
            C15 component of Cij.
        C16 : float
            C16 component of Cij.
        C22 : float
            C22 component of Cij.
        C23 : float
            C23 component of Cij.
        C24 : float
            C24 component of Cij.
        C25 : float
            C25 component of Cij.
        C26 : float
            C26 component of Cij.
        C33 : float
            C33 component of Cij.
        C34 : float
            C34 component of Cij.
        C35 : float
            C35 component of Cij.
        C36 : float
            C36 component of Cij.
        C44 : float
            C44 component of Cij.
        C45 : float
            C45 component of Cij.
        C46 : float
            C46 component of Cij.
        C55 : float
            C55 component of Cij.
        C56 : float
            C56 component of Cij.
        C66 : float
            C66 component of Cij.
        """
        if (C11 is None or C12 is None or C13 is None or
            C14 is None or C15 is None or C16 is None or
            C22 is None or C23 is None or C24 is None or
            C25 is None or C26 is None or C33 is None or
            C34 is None or C35 is None or C36 is None or
            C44 is None or C45 is None or C46 is None or
            C55 is None or C56 is None or C66 is None):
            raise TypeError('triclinic style requires all 6x6 Cij where i <= j')

        # Build Cij array
        self.Cij = np.array([[C11, C12, C13, C14, C15, C16],
                             [C12, C22, C23, C24, C25, C26],
                             [C13, C23, C33, C34, C35, C36],
                             [C14, C24, C34, C44, C45, C46],
                             [C15, C25, C35, C45, C55, C56],
                             [C16, C26, C36, C46, C56, C66]])

    def normalized_as(self, crystal_system: str) -> ElasticConstants:
        """
        Returns a new ElasticConstants object where values of the current are
        averaged or zeroed out according to a standard crystal system setting.
        NOTE: no validation checks are made to evaluate whether such
        normalizations should be done! That is left up to you (compare values
        before and after normalization).
        
        Parameters
        ----------
        crystal_system : str
            Indicates the crystal system representation to use when building a
            data model.
            
        Returns
        -------
        atomman.ElasticConstants
            The elastic constants normalized according to the crystal system
            symmetries.
        """

        c = self.Cij
        c_dict = {}

        if crystal_system == 'isotropic':
            c_dict['mu'] = self.shear()
            c_dict['K'] = self.bulk()

        elif crystal_system == 'cubic':
            c_dict['C11'] = (c[0,0] + c[1,1] + c[2,2]) / 3
            c_dict['C12'] = (c[0,1] + c[0,2] + c[1,2]) / 3
            c_dict['C44'] = (c[3,3] + c[4,4] + c[5,5]) / 3

        elif crystal_system == 'hexagonal':
            c_dict['C11'] = (c[0,0] + c[1,1]) / 2
            c_dict['C33'] = c[2,2]
            c_dict['C12'] = (c[0,1] + (c[0,0] - 2*c[5,5])) / 2
            c_dict['C13'] = (c[0,2] + c[1,2]) / 2
            c_dict['C44'] = (c[3,3] + c[4,4]) / 2

        elif crystal_system == 'tetragonal':
            c_dict['C11'] = (c[0,0] + c[1,1]) / 2
            c_dict['C33'] = c[2,2]
            c_dict['C12'] = c[0,1]
            c_dict['C13'] = (c[0,2] + c[1,2]) / 2
            c_dict['C16'] = (c[0,5] - c[1,5]) / 2
            c_dict['C44'] = (c[3,3] + c[4,4]) / 2
            c_dict['C66'] = c[5,5]

        elif crystal_system == 'rhombohedral':
            c_dict['C11'] = (c[0,0] + c[1,1]) / 2
            c_dict['C33'] = c[2,2]
            c_dict['C12'] = (c[0,1] + (c[0,0] - 2*c[5,5])) / 2
            c_dict['C13'] = (c[0,2] + c[1,2]) / 2
            c_dict['C14'] = (c[0,3] - c[1,3]) / 2
            c_dict['C15'] = (c[0,4] - c[1,4] - c[3,5]) / 3
            c_dict['C44'] = (c[3,3] + c[4,4]) / 2

        elif crystal_system == 'orthorhombic':
            c_dict['C11'] = c[0,0]
            c_dict['C22'] = c[1,1]
            c_dict['C33'] = c[2,2]
            c_dict['C12'] = c[0,1]
            c_dict['C13'] = c[0,2]
            c_dict['C23'] = c[1,2]
            c_dict['C44'] = c[3,3]
            c_dict['C55'] = c[4,4]
            c_dict['C66'] = c[5,5]

        elif crystal_system == 'monoclinic':
            c_dict['C11'] = c[0,0]
            c_dict['C22'] = c[1,1]
            c_dict['C33'] = c[2,2]
            c_dict['C12'] = c[0,1]
            c_dict['C13'] = c[0,2]
            c_dict['C15'] = c[0,4]
            c_dict['C23'] = c[1,2]
            c_dict['C25'] = c[1,4]
            c_dict['C35'] = c[2,4]
            c_dict['C46'] = c[3,5]
            c_dict['C44'] = c[3,3]
            c_dict['C55'] = c[4,4]
            c_dict['C66'] = c[5,5]

        elif crystal_system == 'triclinic':
            c_dict['Cij'] = c

        else:
            raise ValueError('Invalid crystal_system: ' + crystal_system)

        return ElasticConstants(**c_dict)

    def is_normal(self,
                  crystal_system: str,
                  atol: float = 1e-4,
                  rtol: float = 1e-4) -> bool:
        """
        Checks if current elastic constants agree with values normalized to
        a specified crystal family (within tolerances).
        
        Parameters
        ----------
        crystal_system : str
            Indicates the crystal system representation to use when building a
            data model.
        atol : float, optional
            Absolute tolerance to use.  Default value is 1e-4.
        rtol : float, optional
            Relative tolerance to use.  Default value is 1e-4.
        
        Returns
        -------
        bool
            True if all Cij match within the tolerances, false otherwise.
        """
        return np.allclose(self.Cij, self.normalized_as(crystal_system).Cij,
                           atol=atol, rtol=rtol)

    def model(self,
              model: Union[str, io.IOBase, DM, None] = None,
              unit: Optional[str] = None,
              crystal_system: str = 'triclinic') -> Optional[DM]:
        """
        Return or set DataModelDict representation of the elastic constants.
        
        Parameters
        ----------
        model : DataModelDict, string, or file-like object, optional
            Data model containing exactly one 'elastic-constants' branch to
            read.
        unit : str, optional
            Units or pressure to save values in when building a data model.
            Default value is None (no conversion).
        crystal_system : str, optional
            Indicates the crystal system representation to normalize by.
            Default value is 'triclinic', i.e. no normalization.
        
        Returns
        -------
        DataModelDict
            If model is not given as a parameter.
        """

        # Set values if model given
        if model is not None:

            # Find elastic-constants element
            model = DM(model).find('elastic-constants')

            # Read in values
            if 'Cij' in model:
                # New format
                self.Cij = uc.value_unit(model['Cij'])
            else:
                # Old format
                c_dict = {}
                for C in model['C']:
                    key = 'C' + C['ij'][0] + C['ij'][2]
                    c_dict[key] = uc.value_unit(C['stiffness'])
                self.Cij = ElasticConstants(**c_dict).Cij

        # Return DataModelDict if model not given
        else:
            normCij = self.normalized_as(crystal_system).Cij
            model = DM()
            model['elastic-constants'] = DM()
            model['elastic-constants']['Cij'] = uc.model(normCij, unit)

            return model

    def bulk(self, style: str = 'Hill') -> float:
        """
        Returns a bulk modulus estimate.
        
        Parameters
        ----------
        style : str
            Indicates which style of estimate to use.  Default value is 'Hill'.
            - 'Hill' -- Hill estimate (average of Voigt and Reuss).
            - 'Voigt' -- Voigt estimate. Uses Cij.
            - 'Reuss' -- Reuss estimate. Uses Sij.
        """
        if style == 'Hill':
            return (self.bulk('Voigt') + self.bulk('Reuss')) / 2

        elif style == 'Voigt':
            c = self.Cij
            return ( (c[0,0] + c[1,1] + c[2,2]) + 2*(c[0,1] + c[1,2] + c[0,2]) ) / 9

        elif style == 'Reuss':
            s = self.Sij
            return 1 / ( (s[0,0] + s[1,1] + s[2,2]) + 2*(s[0,1] + s[1,2] + s[0,2]) )

        else:
            raise ValueError('Unknown estimate style')

    def shear(self, style: str = 'Hill') -> float:
        """
        Returns a shear modulus estimate.
        
        Parameters
        ----------
        style : str
            Indicates which style of estimate to use.  Default value is 'Hill'.
            - 'Hill' -- Hill estimate (average of Voigt and Reuss).
            - 'Voigt' -- Voigt estimate. Uses Cij.
            - 'Reuss' -- Reuss estimate. Uses Sij.
        """
        if style == 'Hill':
            return (self.shear('Voigt') + self.shear('Reuss')) / 2

        elif style == 'Voigt':
            c = self.Cij
            return ( (c[0,0] + c[1,1] + c[2,2]) - (c[0,1] + c[1,2] + c[0,2]) + 3*(c[3,3] + c[4,4] + c[5,5]) ) / 15

        elif style == 'Reuss':
            s = self.Sij
            return 15 / ( 4*(s[0,0] + s[1,1] + s[2,2]) - 4*(s[0,1] + s[1,2] + s[0,2]) + 3*(s[3,3] + s[4,4] + s[5,5]) )

        else:
            raise ValueError('Unknown estimate style')
