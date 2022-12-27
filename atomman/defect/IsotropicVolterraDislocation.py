# coding: utf-8
# Standard Python libraries
import warnings
from typing import Optional, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
from . import VolterraDislocation
from .. import Box, ElasticConstants

class IsotropicVolterraDislocation(VolterraDislocation):
    """
    Class for representing the isotropic Volterra solution for a straight dislocation.
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
        Computes the elastic solution for an isotropic volterra dislocation.

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
        # Check that C is isotropic
        if not C.is_normal('isotropic', atol=0.0, rtol=1e-4):
            raise ValueError('C must be isotropic elastic constants')
        C = C.normalized_as('isotropic')

        VolterraDislocation.solve(self, C, burgers, ξ_uvw=ξ_uvw, slip_hkl=slip_hkl,
                                  transform=transform, axes=axes, box=box,
                                  m=m, n=n, cart_axes=cart_axes, tol=tol)

        # Save shear modulus and Poissons ratio
        bulk = self.C.bulk()
        self.__mu = self.C.shear()
        self.__nu = (3 * bulk - 2 * self.mu) / (2 * (3 * bulk + self.mu))

    @property
    def mu(self) -> float:
        """float: The isotropic shear modulus"""
        return self.__mu

    @property
    def nu(self) -> float:
        """float: The isotropic Poisson's ratio"""
        return self.__nu

    @property
    def K_tensor(self):
        """numpy.ndarray : The energy coefficient tensor"""

        # Construct K_tensor in standard setting
        K_e = self.mu / (1 - self.nu)
        K_s = self.mu
        K = np.array([[K_e, 0.0, 0.0],
                      [0.0, K_e, 0.0],
                      [0.0, 0.0, K_s]])

        # Transform tensor from m, n, ξ system
        trans = np.array([self.m, self.n, self.ξ])
        K = trans.T.dot(K.dot(trans))

        # Round away near-zero terms
        K[np.isclose(K / K.max(), 0.0, atol=self.tol)] = 0.0
        return K

    def theta(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Computes arctan(y / x) ranging from -π to π.

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The theta angles for each 3D position.
        """
        pos = np.asarray(pos, dtype=float)
        x = pos.dot(self.m)
        y = pos.dot(self.n)

        # Compute arctan(y/x) for all values
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            theta = np.arctan(y / x)

        # Handle special cases and ensure value range -π to π
        theta[(x == 0) & (y > 0)] = np.pi / 2
        theta[(x == 0) & (y < 0)] = -np.pi / 2
        theta[(x < 0)] += np.pi 
        theta[(theta >= np.pi)] -= 2 * np.pi

        return theta

    def displacement(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent isotropic displacements.

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3D vector displacements at all given points.
        """
        pos = np.asarray(pos)
        if pos.shape == (3,):
            pos = pos.reshape(1,3)

        # Split pos, burgers into components
        x = pos.dot(self.m)
        y = pos.dot(self.n)
        b_s = self.burgers.dot(self.ξ)
        b_e = self.burgers.dot(self.m)
        nu = self.nu

        # Compute displacement components in m, n, ξ directions
        disp_m = b_e / (2 * np.pi) * (self.theta(pos) + (x * y) / (2 * (1 - nu) * (x**2 + y**2)))

        disp_n = b_e / (2 * np.pi) * (-(1 - 2 * nu) / (4 * (1 - nu)) * np.log(x**2 + y**2)
                           + (y**2) / (2 * (1 - nu) * (x**2 + y**2)))

        disp_ξ = b_s / (2 * np.pi) * (self.theta(pos))

        # Combine into array
        disp = np.outer(disp_ξ, self.ξ) + np.outer(disp_m, self.m) + np.outer(disp_n, self.n)

        if disp.shape[0] == 1:
            return disp[0]
        else:
            return disp

    def strain(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent isotropic strains.  The equations used are derived
        from ϵ_ij = S_ijkl σ_kl.

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3x3 strain states at all given points.
        """
        pos = np.asarray(pos)
        if pos.shape == (3,):
            pos = pos.reshape(1,3)

        # Split pos, burgers into components
        x = pos.dot(self.m)
        y = pos.dot(self.n)
        b_s = self.burgers.dot(self.ξ)
        b_e = self.burgers.dot(self.m)
        nu = self.nu

        # Initialize empty strain array
        strain = np.empty(pos.shape[:-1] + (3,3))

        # Strain components due to b_s
        strain[..., 0, 2] = strain[..., 2, 0] = -b_s * y / (4 * np.pi * (x**2 + y**2))
        strain[..., 1, 2] = strain[..., 2, 1] =  b_s * x / (4 * np.pi * (x**2 + y**2))

        # Shear strain components due to b_e
        strain[..., 0, 1] = strain[..., 1, 0] = b_e * (x * (x**2 - y**2)) / (4 * np.pi * (1 - nu) * (x**2 + y**2)**2)

        # Normal strain components due to b_e
        common = b_e * y / (4 * np.pi * (1 - nu**2) * (x**2 + y**2)**2)
        strain[..., 0, 0] = common * (x**2 * (-3 - nu + 2 * nu**2) + y**2 * (-1 + nu + 2 * nu**2))
        strain[..., 1, 1] = common * (x**2 * (1 + 3 * nu + 2 * nu**2) + y**2 * (-1 + nu + 2 * nu**2))

        # Set zero values
        strain[..., 2, 2] = 0.0

        # Get the reverse transformation matrix
        transform = np.array([self.m, self.n, self.ξ]).T

        # Transform strains
        strain = np.einsum('mi, nj, ...ij -> ...mn', transform, transform, strain)

        if strain.shape[0] == 1:
            return strain[0]
        else:
            return strain

    def stress(self, pos: npt.ArrayLike) -> np.ndarray:
        """
        Compute the position-dependent isotropic stresses.

        Parameters
        ----------
        pos : array-like object
            3D vector position(s).

        Returns
        -------
        numpy.ndarray
            The computed 3x3 stress states at all given points.
        """
        pos = np.asarray(pos)
        if pos.shape == (3,):
            pos = pos.reshape(1,3)

        # Split pos, burgers into components
        x = pos.dot(self.m)
        y = pos.dot(self.n)
        b_s = self.burgers.dot(self.ξ)
        b_e = self.burgers.dot(self.m)
        nu = self.nu
        mu = self.mu

        # Initialize empty stress array
        stress = np.empty(pos.shape[:-1] + (3,3))

        # Stress components due to b_s
        pre_s = mu * b_s / (2 * np.pi)
        stress[..., 0, 2] = stress[..., 2, 0] =-pre_s * y / (x**2 + y**2)
        stress[..., 1, 2] = stress[..., 2, 1] = pre_s * x / (x**2 + y**2)

        # Stress components due to b_e
        pre_e = mu * b_e / (2 * np.pi * (1 - nu))
        stress[..., 0, 0] =-pre_e * (y * (3 * x**2 + y**2)) / (x**2 + y**2)**2
        stress[..., 1, 1] = pre_e * (y * (x**2 - y**2)) / (x**2 + y**2)**2
        stress[..., 2, 2] = nu * (stress[..., 0, 0] + stress[..., 1, 1])
        stress[..., 0, 1] = stress[..., 1, 0] = pre_e * (x * (x**2 - y**2)) / (x**2 + y**2)**2

        # Get the reverse transformation matrix
        transform = np.array([self.m, self.n, self.ξ]).T

        # Transform stresses
        stress = np.einsum('mi, nj, ...ij -> ...mn', transform, transform, stress)

        if stress.shape[0] == 1:
            return stress[0]
        else:
            return stress
