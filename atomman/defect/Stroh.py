# coding: utf-8
# Standard Python libraries
from copy import deepcopy
from typing import Optional, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from . import VolterraDislocation
from .. import Box, ElasticConstants

class Stroh(VolterraDislocation):
    """
    Class for solving the Eshelby anisotropic solution for a straight
    dislocation or crack using the Stroh method.
    """

    def solve(self,
              C: ElasticConstants,
              burgers: npt.ArrayLike,
              ξ_uvw: Optional[npt.ArrayLike] = None,
              slip_hkl: Optional[npt.ArrayLike] = None,
              transform: Optional[npt.ArrayLike] = None,
              axes: Optional[npt.ArrayLike] = None,
              box: Optional[Box] = None,
              m: Union[str, npt.ArrayLike] = 'x',
              n: Union[str, npt.ArrayLike] = 'y',
              cart_axes: bool = False,
              tol: float = 1e-8):
        """
        Computes the elastic solution for an anisotropic volterra dislocation.

        Parameters
        ----------
        C : atomman.ElasticConstants
            The medium's elastic constants.
        burgers : array-like object
            The dislocation's Burgers vector.
        ξ_uvw : array-like object
            The Miller crystal vector associated with the dislocation's line
            direction.  Must be given with slip_hkl to identify the
            transformation matrix to use on C and burgers.
        slip_hkl : array-like object
            The Miller plane indices associated with the dislocation's slip
            plane.  Must be given with slip_hkl to identify the
            transformation matrix to use on C and burgers.
        transform : array-like object, optional
            A 3x3 set of orthogonal Cartesian vectors that define the
            transformation matrix to use on C and burgers to convert from the
            standard (unit cell) and dislocation orientations.  The 3 vectors
            will automatically be converted into unit vectors.  Using this is
            an alternative to using ξ_uvw and slip_hkl.
        axes : array-like object, optional
            Same as transform.  Retained for backwards compatibility.
        box : atomman.Box, optional
            The unit cell's box that crystal vectors are taken with respect to.
            If not given, will use a cubic box with a=1 meaning that burgers,
            ξ_uvw and slip_hkl will be interpreted as Cartesian vectors.
        m : str or array-like object, optional
            The 3D Cartesian unit vector to align with the dislocation solution's m-axis,
            i.e. the in-plane direction perpendicular to the dislocation line.  Also
            accepts str values of 'x', 'y', or 'z', in which case the dislocation axis will
            be aligned with the corresponding Cartesian axis.  Default value is 'x'.
        n : str or array-like object, optional
            The 3D Cartesian unit vector to align with the dislocation solution's n-axis,
            i.e. the slip plane normal. Also accepts str values of 'x', 'y', or 'z', in
            which case the dislocation axis will be aligned with the corresponding Cartesian
            axis. Default value is 'y'.
        cart_axes : bool, optional
            Setting this to True will also perform an assertion check that the m- and n-axes
            are both aligned with Cartesian axes. This is a requirement for some of the
            atomic configuration generators. Default value is False as the elastic solution
            by itself does not require the limitation.
        tol : float
            Tolerance parameter used to round off near-zero values.  Default
            value is 1e-8.
        """
        VolterraDislocation.solve(self, C, burgers, ξ_uvw=ξ_uvw, slip_hkl=slip_hkl,
                                  transform=transform, axes=axes, box=box,
                                  m=m, n=n, cart_axes=cart_axes, tol=tol)

        # Pull out full 3x3x3x3 elastic constants matrix
        Cijkl = self.C.Cijkl

        # Matrices of Cijkl constants used to construct N
        mm = np.einsum('i,ijkl,l', self.m, Cijkl, self.m)
        mn = np.einsum('i,ijkl,l', self.m, Cijkl, self.n)
        nm = np.einsum('i,ijkl,l', self.n, Cijkl, self.m)
        nn = np.einsum('i,ijkl,l', self.n, Cijkl, self.n)

        # The four 3x3 matrices that represent the quadrants of N
        NB = -np.linalg.inv(nn)
        NA = NB.dot(nm)
        NC = mn.dot(NA) + mm
        ND = mn.dot(NB)

        # N is the 6x6 array, where the eigenvalues are the roots p
        # and the eigenvectors give A and L
        N =  np.array(np.vstack((np.hstack((NA, NB)), np.hstack((NC, ND)))))
        
        # Calculate the eigenvectors and eigenvalues
        eig = np.linalg.eig(N)
        p = eig[0]
        eigvec = np.transpose(eig[1])

        # Separate the eigenvectors into A and L
        A = np.array([eigvec[0,:3], eigvec[1,:3], eigvec[2,:3], eigvec[3,:3], eigvec[4,:3], eigvec[5,:3]])
        L = np.array([eigvec[0,3:], eigvec[1,3:], eigvec[2,3:], eigvec[3,3:], eigvec[4,3:], eigvec[5,3:]])

        # Calculate normalization factor k
        k = 1. / (2. * np.einsum('si,si->s', A, L))

        # Calculation verification checks
        try:
            assert np.allclose(np.einsum('s,si,sj->ij', k,A,L), np.identity(3, dtype='complex128'), atol=tol)
            assert np.allclose(np.einsum('s,si,sj->ij', k,A,A), np.zeros((3,3), dtype='complex128'), atol=tol)
            assert np.allclose(np.einsum('s,si,sj->ij', k,L,L), np.zeros((3,3), dtype='complex128'), atol=tol)
            assert np.allclose(np.einsum('s,t,si,ti->st', k**.5,k**.5,A,L) + np.einsum('s,t,ti,si->st', k**.5,k**.5,A,L),
                                         np.identity(6, dtype='complex128'), atol = tol)
        except AssertionError as err:
            raise ValueError('Stroh checks failed!') from err

        # Assign property values
        self.__p = p
        self.__A = A
        self.__L = L
        self.__k = k

        # Check that K_tensor is real
        if self.K_tensor.dtype == 'complex128':
            raise ValueError('Solution not real: check elastic constants')

    @property
    def p(self) -> np.ndarray:
        """numpy.ndarray : p eigenvalues"""
        return deepcopy(self.__p)

    @property
    def A(self) -> np.ndarray:
        """numpy.ndarray : A eigenvectors"""
        return deepcopy(self.__A)

    @property
    def L(self) -> np.ndarray:
        """numpy.ndarray : L eigenvectors"""
        return deepcopy(self.__L)

    @property
    def k(self) -> np.ndarray:
        """numpy.ndarray : k normalization factors"""
        return deepcopy(self.__k)

    @property
    def K_tensor(self) -> np.ndarray:
        """numpy.ndarray : The energy coefficient tensor"""

        # ii is imaginary unit
        ii = np.array([1.j])

        # updn is alternating +-
        updn = np.array([1, -1, 1, -1, 1, -1])

        # Compute K_tensor
        K = ii * np.einsum('s,s,si,sj->ij', updn, self.k, self.L, self.L)

        # Round away near-zero terms
        K = np.real_if_close(K, tol=self.tol)
        K[np.isclose(K / K.max(), 0.0, atol=self.tol)] = 0.0

        return K

    def displacement(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent anisotropic displacement.

            u_i = 1 / (2 π i) (Σ_a +- k_a A_ai (L_aj*burgers_j) ln(η_a))

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3D vector displacements at all given points.
        """
        # ii is imaginary unit
        ii = np.array([1.j])

        # updn is alternating +-
        updn = np.array([1, -1, 1, -1, 1, -1])

        # Compute dot of L and b
        Lb = self.L.dot(self.burgers)

        # Combine k, updn, Lb
        kLb = self.k * updn * Lb

        # Convert pos to eta coordinate
        eta = self.eta(pos)

        # Compute the displacements
        disp = 1 / (2 * np.pi * ii) * np.einsum('a, ai, ...a -> ...i',
                                                kLb, self.A, np.log(eta))

        # Round away near-zero terms
        disp = np.real_if_close(disp, tol=self.tol)

        # Reduce single-value solutions
        if disp.shape[0] == 1:
            return disp[0]
        else:
            return disp

    def strain(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent anisotropic stresses.

            ϵ_ij = 1 / (4 π i) (Σ_a +- k_a (mpn_ai A_aj + mpn_aj A_ai) (L_am * burgers_m) / η_a)

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3x3 stress states at all given points.
        """
        # ii is imaginary unit
        ii = np.array([1.j])

        # updn is alternating +-
        updn = np.array([1, -1, 1, -1, 1, -1])

        # Compute dot of L and b
        Lb = self.L.dot(self.burgers)

        # Combine k, updn, Lb
        kLb = self.k * updn * Lb

        # Compute symmetric Ampn factor
        mpn = self.m + np.outer(self.p, self.n)
        Ampn = np.einsum('ai, aj -> aij', mpn, self.A) + np.einsum('aj, ai -> aij', mpn, self.A)

        # Convert pos to eta coordinate
        eta = self.eta(pos)

        # Compute the strains
        strain = 1 / (4 * np.pi * ii) * np.einsum('a, aij, ...a -> ...ij',
                                                  kLb, Ampn, 1/eta)

        # Round away near-zero terms
        strain = np.real_if_close(strain, tol=self.tol)

        if strain.shape[0] == 1:
            return strain[0]
        else:
            return strain

    def stress(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent anisotropic stresses.

            σ_ij = 1 / (2 π i) (Σ_a +- k_a C_ijkl mpn_al A_ak (L_am*burgers_m) / η_a)

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3x3 stress states at all given points.
        """
        # ii is imaginary unit
        ii = np.array([1.j])

        # updn is alternating +-
        updn = np.array([1, -1, 1, -1, 1, -1])

        # Compute dot of L and b
        Lb = self.L.dot(self.burgers)

        # Combine k, updn, Lb
        kLb = self.k * updn * Lb

        # Compute mpn factor
        mpn = self.m + np.outer(self.p, self.n)
        Ampn = np.einsum('ai, aj -> aij', mpn, self.A)

        # Convert pos to eta coordinate
        eta = self.eta(pos)

        # Compute the stresses
        stress = 1 / (2 * np.pi * ii) * np.einsum('a, ijkl, alk, ...a -> ...ij',
                                                  kLb, self.C.Cijkl, Ampn, 1/eta)

        # Round away near-zero terms
        stress = np.real_if_close(stress, tol=self.tol)

        if stress.shape[0] == 1:
            return stress[0]
        else:
            return stress

    def eta(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the eta coordinates based on positions, p, m and n.  Used by
        displacement() and stress().

            η_a = x_i m_i + p_a x_j n_j

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed eta factor at all given points.
        """

        x = np.dot(pos, self.m)
        y = np.dot(pos, self.n)

        return (x + np.outer(self.p, y)).T
