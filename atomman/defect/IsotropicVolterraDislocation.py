# coding: utf-8
# Standard Python libraries
import warnings

# http://www.numpy.org/
import numpy as np

# atomman imports
from . import VolterraDislocation
from ..tools import vect_angle, axes_check

class IsotropicVolterraDislocation(VolterraDislocation):
    """
    Class for representing the isotropic Volterra solution for a straight dislocation.
    """
    
    def solve(self, C, burgers, axes=None, m=[1,0,0], n=[0,1,0], tol=1e-8):
        """
        Sets the parameters for the solution to use.
        
        Parameters
        ----------
        C : atomman.ElasticConstants
            The medium's elastic constants.  Must be isotropic.
        burgers : array-like object
            The dislocation's Cartesian Burgers vector.
        axes : array-like object, optional
            3x3 set of rotational axes for the system. If given, C and burgers
            will be transformed using axes.
        m : array-like object, optional
            The m unit vector for the solution.  m, n, and u (dislocation
            line) should be right-hand orthogonal.  Default value is [1,0,0]
            (x-axis).
        n : array-like object, optional
            The n unit vector for the solution.  m, n, and u (dislocation
            line) should be right-hand orthogonal.  Default value is [0,1,0]
            (y-axis). n is normal to the dislocation slip plane.
        tol : float
            Tolerance parameter used to check for compatibility of the other
            parameters.  Default value is 1e-8.
        """
        VolterraDislocation.solve(self, C, burgers, axes=axes, m=m, n=n, tol=tol)
        
        # Check that C is isotropic
        if not C.is_normal('isotropic', atol=0.0, rtol=1e-4):
            raise ValueError('C must be isotropic elastic constants')
        self.__C = C.normalized_as('isotropic')
        
        # Save shear modulus and Poissons ratio
        bulk = self.C.bulk()
        self.__mu = self.C.shear()
        self.__nu = (3 * bulk - 2 * self.mu) / (2 * (3 * bulk + self.mu))
    
    @property
    def mu(self):
        return self.__mu
    
    @property
    def nu(self):
        return self.__nu
    
    @property
    def K_tensor(self):
        
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

    def theta(self, pos):
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
    
    def displacement(self, pos):
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

    def stress(self, pos):
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
        
        # Define screw / edge prefactors
        pre_s = mu * b_s / (2 * np.pi)
        pre_e = mu * b_e / (2 * np.pi * (1 - nu))

        # Initialize empty stress array
        sigma = np.empty(pos.shape[:-1] + (3,3))
  
        # Stress components due to b_s
        sigma[..., 0, 2] = sigma[..., 2, 0] =-pre_s * y / (x**2 + y**2)
        sigma[..., 1, 2] = sigma[..., 2, 1] = pre_s * x / (x**2 + y**2)

        # Stress components due to b_e
        sigma[..., 0, 0] =-pre_e * (y * (3 * x**2 + y**2)) / (x**2 + y**2)**2
        sigma[..., 1, 1] = pre_e * (y * (x**2 - y**2)) / (x**2 + y**2)**2
        sigma[..., 2, 2] = nu * (sigma[..., 0, 0] + sigma[..., 1, 1])
        sigma[..., 0, 1] = sigma[..., 1, 0] = pre_e * (x * (x**2 - y**2)) / (x**2 + y**2)**2

        # Get the reverse transformation matrix
        transform = np.array([self.m, self.n, self.ξ]).T

        # Transform stresses
        return np.einsum('mi,nj,...ij->...mn', transform, transform, sigma)
