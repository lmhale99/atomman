# coding: utf-8

# Standard Python libraries
import io
import warnings
from typing import Optional, Tuple, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# https://www.scipy.org/
from scipy.optimize import minimize, OptimizeResult

# https://matplotlib.org/
import matplotlib.pyplot as plt

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
import atomman.unitconvert as uc
from . import GammaSurface, VolterraDislocation

class SDVPN(object):
    """
    Class representation of the semidiscrete variational Peierls-Nabarro
    dislocation model.
    """
    
    def __init__(self,
                 volterra: Optional[VolterraDislocation] = None,
                 gamma: Optional[GammaSurface] = None,
                 model: Union[str, io.IOBase, DM, None] = None, 
                 tau: npt.ArrayLike = np.zeros((3,3)),
                 alpha: float = 0.0,
                 beta: npt.ArrayLike = np.zeros((3,3)),
                 cutofflongrange: Optional[float] = None,
                 fullstress: bool = True,
                 cdiffelastic: bool = False,
                 cdiffsurface: bool = True,
                 cdiffstress: bool = False,
                 min_method: str = 'Powell',
                 min_kwargs: Optional[dict] = None,
                 min_options: Optional[dict] = None):
        """
        Initializes an SDVPN object.
        
        Parameters
        ----------
        volterra : atomman.defect.VolterraDislocation, optional
            The elastic solution for a Volterra dislocation to use as the basis
            of the model. Either volterra or model are required, and both cannot
            be given at the same time.
        gamma : atomman.defect.GammaSurface, optional
            The gamma surface to use for the solution.  Required unless model
            is given and the model content contains gamma surface data.
        model : str or DataModelDict, optional
            Saved data from previous SDVPN runs to load.  Either volterra or
            model are required, and both cannot be given at the same time.
        tau : numpy.ndarray, optional
            A (3,3) array giving the stress tensor to apply to the system
            using the stress energy term.  Only the xy, yy, and yz components
            are used.  Default value is all zeros.
        alpha : list of float, optional
            The alpha coefficient(s) used by the nonlocal energy term.  Default
            value is [0.0].
        beta : numpy.ndarray, optional
            The (3,3) array of beta coefficient(s) used by the surface energy
            term.  Default value is all zeros.
        cutofflongrange : float, optional
            The cutoff distance to use for computing the long-range energy.
            Default value is 1000 angstroms.
        fullstress : bool, optional
            Flag indicating which stress energy algorithm to use.  Default
            value is True.
        cdiffelastic : bool, optional
            Flag indicating if the dislocation density for the elastic energy
            component is computed with central difference (True) or simply
            neighboring values (False).  Default value is False.
        cdiffsurface : bool, optional
            Flag indicating if the dislocation density for the surface energy
            component is computed with central difference (True) or simply
            neighboring values (False).  Default value is True.
        cdiffstress : bool, optional
            Flag indicating if the dislocation density for the stress energy
            component is computed with central difference (True) or simply
            neighboring values (False).  Only matters if fullstress is True.
            Default value is False.
        min_method : str, optional
            The scipy.optimize.minimize method to use.  Default value is
            'Powell'.
        min_kwargs : dict, optional
            Any keyword arguments to pass on to scipy.optimize.minimize besides
            the coordinates, method and options.
        min_options : dict, optional
            Any options to pass on to scipy.optimize.minimize. Default value
            is {}.
        """
        
        # Load solution from existing model
        if model is not None:
            if volterra is not None:
                raise ValueError('model cannot be given with volterra')
            
            self.load(model, gamma=gamma)
        
        # Extract parameters and check solution compatibility
        elif volterra is not None:
            
            # Check that gamma is given
            if gamma is None:
                raise ValueError('gamma is required if volterra is given')
            
            # Check that burgers is in the slip plane
            if not np.isclose(np.dot(volterra.n, volterra.burgers), 0.0):
                raise ValueError('dislocation burgers vector must be in the slip plane')
            
            # Get m, n, ξ, K_tensor, burgers, and transform from volterra
            m, n, ξ = volterra.m, volterra.n, volterra.ξ
            K_tensor = volterra.K_tensor
            burgers = volterra.burgers
            transform = volterra.transform
            
            # Transform K_tensor, burgers and transform to [m,n,ξ] orientation
            mnξ = np.array([m, n, ξ]) # This is transformation matrix to [m,n,ξ] setting
            K_tensor = mnξ.dot(K_tensor.dot(mnξ.T))
            burgers = mnξ.dot(burgers)
            transform = np.matmul(mnξ, transform)
            
            # Check if dislocation system is compatible with gamma surface
            planenormal = transform.dot(gamma.planenormal)
            if not np.isclose(planenormal[0], 0.0) or not np.isclose(planenormal[2], 0.0):
                raise ValueError('different slip planes for gamma and volterra found')

            # Set basic solution definition
            self.__K_tensor = K_tensor
            self.__burgers = burgers
            self.__transform = transform
            self.__gamma = gamma

            # Set options
            self.tau = tau
            self.alpha = alpha
            self.beta = beta
            if cutofflongrange is None:
                self.cutofflongrange = uc.set_in_units(1000, 'angstrom')
            else:
                self.cutofflongrange = cutofflongrange
            self.fullstress = fullstress
            self.cdiffelastic = cdiffelastic
            self.cdiffsurface = cdiffsurface
            self.cdiffstress = cdiffstress
            self.min_method = min_method
            self.min_options = min_options
            self.min_kwargs = min_kwargs

        else:
            raise ValueError('either dislsol or model must be given')
    
    @property
    def x(self) -> np.ndarray:
        """numpy.ndarray : The x coordinates."""
        try:
            return self.__x
        except:
            raise AttributeError('x values not set yet')
    
    @x.setter
    def x(self, value: npt.ArrayLike):
        value = np.asarray(value, dtype=float)
        assert value.ndim == 1
        diff = value[1:] - value[:-1]
        assert np.allclose(diff[0], diff), 'x values must be evenly spaced'
        assert diff[0] > 0, 'x values must be in increasing order'
        self.__x = value

    @property
    def disregistry(self) -> np.ndarray:
        """numpy.ndarray : The disregistry vector for each x coordinate."""
        try:
            return self.__disregistry
        except:
            raise AttributeError('disregistry values not set yet')

    @disregistry.setter
    def disregistry(self, value: npt.ArrayLike):
        value = np.asarray(value)
        assert value.ndim == 2 and value.shape[1] == 3, 'invalid disregistry dimensions'
        assert np.allclose(value[:,1], 0.0), 'y (i.e. out-of-plane) component of disregistry not supported'
        self.__disregistry = value

    @property
    def K_tensor(self) -> np.ndarray:
        """numpy.ndarray : Dislocation energy coefficient tensor."""
        return self.__K_tensor

    @property
    def burgers(self) -> np.ndarray:
        """numpy.ndarray : Burgers vector."""
        return self.__burgers
        
    @property
    def transform(self) -> np.ndarray:
        """numpy.ndarray : Transformation matrix from standard crystal setting to dislocation solution setting."""
        return self.__transform
    
    @property
    def gamma(self) -> GammaSurface:
        """atomman.defect.GammaSurface : The stacking fault map."""
        return self.__gamma
    
    @property
    def tau(self) -> np.ndarray:
        """numpy.ndarray : The applied 3x3 stress tensor."""
        return self.__tau

    @tau.setter
    def tau(self, value: npt.ArrayLike):
        value = np.asarray(value, dtype=float)
        assert value.shape == (3,3)
        self.__tau = value
    
    @property
    def alpha(self) -> tuple:
        """tuple of float : Coefficients for nonlocal energy correction."""
        return self.__alpha
    
    @alpha.setter
    def alpha(self, value: Union[float, tuple]):
        try:
            self.__alpha = tuple(value)
        except:
            self.__alpha = (value,)

    @property
    def beta(self) -> np.ndarray:
        """numpy.ndarray : 3x3 coefficients for gradient energy correction."""
        return self.__beta
    
    @beta.setter
    def beta(self, value: npt.ArrayLike):
        value = np.asarray(value, dtype=float)
        assert value.shape == (3,3)
        self.__beta = value

    @property
    def cutofflongrange(self) -> float:
        """float : Cutoff distance for long-range elastic energy."""
        return self.__cutofflongrange
    
    @cutofflongrange.setter
    def cutofflongrange(self, value: float):
        self.__cutofflongrange = float(value)

    @property
    def fullstress(self) -> bool:
        """bool : Flag indicating which stress algorithm was used."""
        return self.__fullstress

    @fullstress.setter
    def fullstress(self, value: bool):
        assert isinstance(value, bool)
        self.__fullstress = value
    
    @property
    def cdiffelastic(self) -> bool:
        """bool : Flag indicating if elastic energy used central difference for computing the dislocation density."""
        return self.__cdiffelastic
    
    @cdiffelastic.setter
    def cdiffelastic(self, value: bool):
        assert isinstance(value, bool)
        self.__cdiffelastic = value

    @property
    def cdiffsurface(self) -> bool:
        """bool : Flag indicating if surface energy used central difference for computing the dislocation density."""
        return self.__cdiffsurface
    
    @cdiffsurface.setter
    def cdiffsurface(self, value: bool):
        assert isinstance(value, bool)
        self.__cdiffsurface = value

    @property
    def cdiffstress(self) -> bool:
        """bool : Flag indicating if stress energy used central difference for computing the dislocation density."""
        return self.__cdiffstress
    
    @cdiffstress.setter
    def cdiffstress(self, value: bool):
        assert isinstance(value, bool)
        self.__cdiffstress = value

    @property
    def min_method(self) -> str:
        """str : scipy.optimize.minimize method used."""
        return self.__min_method
    
    @min_method.setter
    def min_method(self, value: str):
        self.__min_method = str(value)

    @property
    def min_options(self) -> dict:
        """dict : scipy.optimize.minimize options used."""
        return self.__min_options

    @min_options.setter
    def min_options(self, value: Optional[dict]):
        if value is None:
            self.__min_options = {}
        elif isinstance(value, dict):
            self.__min_options = value
        else:
            raise TypeError('min_options must be a dict')

    @property
    def min_kwargs(self) -> dict:
        """dict : scipy.optimize.minimize keywords used."""
        return self.__min_kwargs
    
    @min_kwargs.setter
    def min_kwargs(self, value: Optional[dict]):
        if value is None:
            self.__min_kwargs = {}
        elif isinstance(value, dict):
            self.__min_kwargs = value
        else:
            raise TypeError('min_kwargs must be a dict')

    @property
    def res(self) -> OptimizeResult:
        """OptimizeResult : scipy.optimize.minimize result."""
        try:
            return self.__res
        except:
            return None
    
    def solve(self,
              x: Optional[npt.ArrayLike] = None,
              disregistry: Optional[npt.ArrayLike] = None,
              tau: Optional[npt.ArrayLike] = None,
              alpha: Optional[list] = None,
              beta: Optional[npt.ArrayLike] = None,
              cutofflongrange: Optional[float] = None,
              fullstress: Optional[bool] = None,
              cdiffelastic: Optional[bool] = None,
              cdiffsurface: Optional[bool] = None,
              cdiffstress: Optional[bool] = None,
              min_method: Optional[str] = None,
              min_kwargs: Optional[dict] = None,
              min_options: Optional[dict] = None):
        """
        Solves the semidiscrete variational Peierls-Nabarro dislocation
        disregistry through energy minimization using the set class
        properties.  All parameters are optional keyword arguments that
        can be used to change any of the previous settings.
        
        Parameters
        ----------
        x : numpy.ndarray, optional
            An array of shape (N) giving the x coordinates corresponding to
            the disregistry solution.
        disregistry : numpy.ndarray, optional
            A (N,3) array giving the initial disregistry vector guess at each
            x coordinate.
        tau : numpy.ndarray, optional
            A (3,3) array giving the stress tensor to apply to the system
            using the stress energy term.  Only the xy, yy, and yz components
            are used.
        alpha : list of float, optional
            The alpha coefficient(s) used by the nonlocal energy term.
        beta : numpy.ndarray, optional
            The (3,3) array of beta coefficient(s) used by the surface energy
            term.
        cutofflongrange : float, optional
            The cutoff distance to use for computing the long-range energy.
        fullstress : bool, optional
            Flag indicating which stress energy algorithm to use.
        cdiffelastic : bool, optional
            Flag indicating if the dislocation density for the elastic energy
            component is computed with central difference (True) or simply
            neighboring values (False).
        cdiffsurface : bool, optional
            Flag indicating if the dislocation density for the surface energy
            component is computed with central difference (True) or simply
            neighboring values (False).
        cdiffstress : bool, optional
            Flag indicating if the dislocation density for the stress energy
            component is computed with central difference (True) or simply
            neighboring values (False).  Only matters if fullstress is True.
        min_method : str, optional
            The scipy.optimize.minimize method to use.
        min_kwargs : dict, optional
            Any keyword arguments to pass on to scipy.optimize.minimize besides
            the coordinates, method and options.
        min_options : dict, optional
            Any options to pass on to scipy.optimize.minimize.
        """

        # Change attribute values if given
        if x is not None:
            self.x = x
        if disregistry is not None:
            self.disregistry = disregistry
        if tau is not None:
            self.tau = tau
        if alpha is not None:
            self.alpha = alpha
        if beta is not None:
            self.beta = beta
        if cutofflongrange is not None:
            self.cutofflongrange = cutofflongrange
        if fullstress is not None:
            self.fullstress = fullstress
        if cdiffelastic is not None:
            self.cdiffelastic = cdiffelastic
        if cdiffsurface is not None:
            self.cdiffsurface = cdiffsurface
        if cdiffstress is not None:
            self.cdiffstress = cdiffstress
        if min_method is not None:
            self.min_method = min_method
        if min_options is not None:
            self.min_options = min_options
        if min_kwargs is not None:
            self.min_kwargs = min_kwargs

        # Check that x and disregistry exist and are of the same length
        if len(self.x) != len(self.disregistry):
            raise ValueError('x and disregistry are not of the same length')
        
        # Define subfunctions
        def decompose(d):
            """Breaks disregistry into components for minimization"""
            d13 = np.concatenate([d[1:-1, 0], d[1:-1, 2]])
            first = d[0]
            last = d[-1]
            return d13, first, last

        def recompose(d13, first, last):
            """Reassembles the disregistry components"""
            half = int(len(d13)/2)
            d = np.zeros((half+2, 3))
            d[0] = first
            d[-1] = last
            d[1:-1, 0] = d13[:half]
            d[1:-1, 2] = d13[half:]
            return d

        def min_func(d13, first, last):
            """Function for minimizing"""
            disregistry = recompose(d13, first, last)
            return self.total_energy(disregistry=disregistry)
        
        # Solve disregistry
        d13, first, last = decompose(self.disregistry)
        res = minimize(min_func, d13, args=(first, last),
                    method=self.min_method, options=self.min_options, **self.min_kwargs)
        self.disregistry = recompose(res.x, first, last)

        self.__res = res
    
    def disldensity(self,
                    x: Optional[npt.ArrayLike] = None,
                    disregistry: Optional[npt.ArrayLike] = None,
                    cdiff: bool = False
                    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes the dislocation density as the numerical derivative of
        disregistry with respect to x.  Uses either neighboring values
        
            ρ[i] = (δ[i] - δ[i-1]) / (x[i] - x[i-1])
        
        or central difference
        
            ρ[i] = (δ[i+1] - δ[i-1]) / (x[i+1] - x[i-1])
        
        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        cdiff : bool, optional
            Flag indicating how to compute the derivative.  A value of False
            (default) compares nearest values while a value of True compares
            next-nearest (i.e., central difference).
        
        Returns
        -------
        newx : numpy.array
            The x positions corresponding to the dislocation density values.
        rho : numpy.array
            The computed dislocation density.
        """
        # Default values are class properties
        if x is None:
            x = self.x
        if disregistry is None:
            disregistry = self.disregistry
        
        # Extract values
        δ = disregistry
        
        if cdiff is False:
            # ρ[i] = (δ[i] - δ[i-1]) / (x[i] - x[i-1])
            ρ = ((δ[1:] - δ[:-1]).T / (x[1:] - x[:-1])).T
            
            # newx is all x except the first
            newx = x[1:]
        
        elif cdiff is True:
            # ρ[i] = (δ[i+1] - δ[i-1]) / (x[i+1] - x[i-1])
            ρ = ((δ[2:] - δ[:-2]).T / (x[2:] - x[:-2])).T
            
            # newx is all x except the first and last
            newx = x[1:-1]
        else:
            raise TypeError('cdiff must be bool')
            
        return (newx, ρ)
    
    def misfit_energy(self,
                      x: Optional[npt.ArrayLike] = None,
                      disregistry: Optional[npt.ArrayLike] = None) -> float:
        """
        Computes the misfit energy for the disregistry using the stored gamma
        surface
        
            E_misfit = Σ γ(δ)Δx
        
        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
            
        Returns
        -------
        float
            The misfit energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.x
        if disregistry is None:
            disregistry = self.disregistry
        
        # Extract values
        δ = disregistry
        Δx = x[1] - x[0]
        transform = self.transform
        gamma = self.gamma
        
        # Strip out y-component of disregistry and transform for gamma
        disreg = np.vstack([δ[:,0], np.zeros(len(δ)), δ[:,2]]).T
        pos = np.inner(disreg, transform.T)
        
        # Σ γ(δ)Δx
        return Δx * gamma.E_gsf(pos=pos).sum()
    
    def elastic_energy(self,
                       x: Optional[npt.ArrayLike] = None,
                       disregistry: Optional[npt.ArrayLike] = None) -> float:
        r"""
        Computes the short-range configuration-dependent elastic energy term
        for the dislocation based on the dislocation density and K_tensor.
        
            E_elastic = 1/(4π) Σ_i Σ_j χ(i,j,Δx) K_lm ρ_l[i] ρ_m[j]
            
            χ(i,j,Δx) = (3/2) Δx² + ψ(i-1,j-1,Δx) + ψ(i,j,Δx) - ψ(i,j-1,Δx) - ψ(j,i-1,Δx)
            
            ψ(i,j,Δx) = (1/2) (i-j)² Δx² ln(\|i-j\|Δx)
        
        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The elastic energy for the dislocation.
        """
        # Define subfunctions
        def χ(i, j, Δx):
            """
            Computes the chi subfunction:
                χ(i,j,Δx) = (3/2) Δx² + ψ(i-1,j-1,Δx) + ψ(i,j,Δx)
                                      - ψ(i,j-1,Δx) - ψ(j,i-1,Δx)
            """
            return 3./2. * Δx**2 + (ψ(i-1, j-1, Δx) + ψ(i, j, Δx)
                                    - ψ(i, j-1, Δx) - ψ(j, i-1, Δx))

        def ψ(i, j, Δx):
            """
            Computes the psi subfunction:
                ψ(i,j,Δx) = (1/2) (i-j)² Δx² ln(|i-j|Δx)
            """
            # Suppress NaN runtime warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                
                # ψ(i,j,Δx) = (1/2) (i-j)² Δx² ln(|i-j|Δx)
                p = 0.5 * (i - j)**2 * Δx**2 * np.log(np.abs(i - j) * Δx)
            
            # Replace NaN values with 0.0
            p[np.isnan(p)] = 0.0
            
            return p
        
        # Default values are class properties
        if x is None:
            x = self.x
        if disregistry is None:
            disregistry = self.disregistry
        
        # Extract values
        δ = disregistry
        Δx = x[1] - x[0]
        cdiff = self.cdiffelastic
        Kij = self.K_tensor
        
        ρ = self.disldensity(x=x, disregistry=δ, cdiff=cdiff)[1]
        j = np.arange(len(ρ), dtype=int)
        
        energy = 0.0
        
        # Compute elastic energy (looping over i, vectorization over j)
        # 1/(4π) Σ_i Σ_j χ(i,j,Δx) K_lm ρ_l[i] ρ_m[j]
        for i in j:
            energy += np.sum( χ(i, j, Δx) * np.inner(ρ[i].dot(Kij), ρ) ) / (4 * np.pi)
        
        return energy
    
    def longrange_energy(self) -> float:
        """
        Computes the long-range elastic energy term for the dislocation using
        the K_tensor, Burgers vector and long-range cutoff.  This term is
        configuration-independent and thus the method takes no parameters.
        
            E_longrange = 1/(2π) K_lm b_l b_m ln(L)
        
        Returns
        -------
        float
            The long-range energy for the dislocation.
        """
        # Extract values
        Kij = self.K_tensor
        b = self.burgers
        L = self.cutofflongrange
        
        # Compute long-range energy
        # 1/(2π) K_lm b_l b_m ln(L)
        return np.inner(b.dot(Kij), b) * np.log(L) / (2 * np.pi)
    
    def stress_energy(self, 
                      x: Optional[npt.ArrayLike] = None,
                      disregistry: Optional[npt.ArrayLike] = None) -> float:
        """
        Computes the stress energy due to the applied stress, tau.
        If fullstress is True, the original stress expression by
        Bulatov and Kaxiras will be used:
        
            E_stress = -1/2 Σ_i (x[i]² - x[i-1]²) ρ_l τ_2l
        
        If fullstress is False, the alternate stress expression by
        Shen and Cheng 10.1016/j.scriptamat.2009.04.047 will be used:
        
            E_stress = -1/2 Σ_i τ_2l (δ_l[i] + δ_l[i+1]) Δx
        
        Note that the Shen and Cheng expression will have a constant
        error associated with it giving an incorrect overall energy, but
        should apply a similar force on the dislocation.

        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The stress energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.x
        if disregistry is None:
            disregistry = self.disregistry
        
        # Extract values
        δ = disregistry
        Δx = x[1] - x[0]
        τ = self.tau
        full = self.fullstress
        cdiff = self.cdiffstress
        
        if full is True:
            ρ = self.disldensity(x=x, disregistry=δ, cdiff=cdiff)[1]
            # -1/2 Σ_i (x[i]² - x[i-1]²) ρ_l τ_2l
            return -0.5 * np.sum((x[1:]**2 - x[:-1]**2) * np.inner(ρ, τ[1,:]))
            
        else:
            # Flip sign on tau so energies match full=True
            τ = -τ
            
            # -1/2 Σ_i τ_2l (δ_l[i] + δ_l[i+1]) Δx
            return -0.5 * np.sum(np.inner(τ[1,:], (δ[:-1] + δ[1:]) * Δx))
    
    def surface_energy(self, 
                       x: Optional[npt.ArrayLike] = None,
                       disregistry: Optional[npt.ArrayLike] = None) -> float:
        """
        Computes the gradient surface energy correction using beta
        coefficients.
        
            E_surface = Σ_j β_lj / 4 Σ_i ρ_l[i]² Δx
        
        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The surface energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.x
        if disregistry is None:
            disregistry = self.disregistry
            
        # Extract values
        δ = disregistry
        Δx = x[1] - x[0]
        β = self.beta
        cdiff = self.cdiffsurface
        
        ρ = self.disldensity(x=x, disregistry=δ, cdiff=cdiff)[1]
        
        # Σ_j β_lj / 4 Σ_i ρ_l[i]² Δx
        return np.sum( np.inner(ρ**2 * Δx, β) ) / 4
    
    def nonlocal_energy(self, 
                        x: Optional[npt.ArrayLike] = None,
                        disregistry: Optional[npt.ArrayLike] = None) -> float:
        """
        Computes the nonlocal energy correction using alpha coefficient(s).
        
            E_nonlocal = Σ_m α_m Σ_i δ[i] (δ[i] - (δ[i+m] + δ[i-m]) / 2) Δx
        
        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The nonlocal energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.x
        if disregistry is None:
            disregistry = self.disregistry
            
        # Extract values
        δ = disregistry
        Δx = x[1] - x[0]
        αs = self.alpha
        
        energy = 0.0
        
        # Σ_m α_m Σ_i δ[i] (δ[i] - (δ[i+m] + δ[i-m]) / 2) Δx
        for num, α in enumerate(αs):
            m = num + 1
            dd = δ[m:-m] - 0.5 * (δ[2*m:] + δ[:-2*m])
            energy += α * np.sum(δ[m:-m] * dd * Δx)
            
        return energy
    
    def total_energy(self, 
                     x: Optional[npt.ArrayLike] = None,
                     disregistry: Optional[npt.ArrayLike] = None) -> float:
        """
        Computes the total energy for the dislocation.
        
            E_total = E_elastic + E_misfit + E_stress + E_surface + E_nonlocal
        
        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        
        Returns
        -------
        float
            The total energy for the dislocation.
        """
        # Default values are class properties
        if x is None:
            x = self.x
        if disregistry is None:
            disregistry = self.disregistry
            
        return (self.misfit_energy(x, disregistry)
                + self.elastic_energy(x, disregistry)
                + self.longrange_energy()
                + self.stress_energy(x, disregistry)
                + self.nonlocal_energy(x, disregistry)
                + self.surface_energy(x, disregistry))
    
    def check_energies(self,
                       x: Optional[npt.ArrayLike] = None,
                       disregistry: Optional[npt.ArrayLike] = None,
                       energyperlength_unit: str = 'eV/Å'):
        """
        Prints a summary string of all computed energy components.

        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        energyperlength_unit : str, optional
            The units of energy per length to report the dislocation line
            energies in.  Default value is 'eV/Å'.
        """
        hasvals = True
        # Default values are class properties
        if x is None:
            try:
                x = self.x
            except:
                hasvals = False
        if disregistry is None:
            try:
                disregistry = self.disregistry
            except:
                hasvals = False

        if hasvals:
            print(f'Dislocation energy terms in {energyperlength_unit}:')
            print('Misfit energy =    ', uc.get_in_units(self.misfit_energy(x, disregistry), energyperlength_unit))
            print('Elastic energy =   ', uc.get_in_units(self.elastic_energy(x, disregistry), energyperlength_unit))
            print('Long-range energy =', uc.get_in_units(self.longrange_energy(), energyperlength_unit))
            print('Stress energy =    ', uc.get_in_units(self.stress_energy(x, disregistry), energyperlength_unit))
            print('Surface energy =   ', uc.get_in_units(self.surface_energy(x, disregistry), energyperlength_unit))
            print('Nonlocal energy =  ', uc.get_in_units(self.nonlocal_energy(x, disregistry), energyperlength_unit))
            print('Total energy =     ', uc.get_in_units(self.total_energy(x, disregistry), energyperlength_unit))
        else:
            print('x and disregistry must be set/given to check energies')
    
    def disregistry_plot(self, 
                         x: Optional[npt.ArrayLike] = None,
                         disregistry: Optional[npt.ArrayLike] = None,
                         figsize: Optional[tuple] = None,
                         length_unit: str = 'Å') -> plt.figure:
        """
        Creates a simple matplotlib figure showing the disregistry profiles.

        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.
        figsize : tuple, optional
            matplotlib figure figsize parameter.  Default value is (10, 6).
        length_unit : str, optional
            The unit of length to display positions and disregistry values in.
            Default value is 'Å'.
            
        Returns
        -------
        matplotlib.pyplot.figure
            The generated figure allowing users to perform additional
            modifications.
        """
        hasvals = True
        # Default values are class properties
        if x is None:
            try:
                x = self.x
            except:
                hasvals = False
        if disregistry is None:
            try:
                disregistry = self.disregistry
            except:
                hasvals = False
        if hasvals:

            if figsize is None:
                figsize = (10, 6)

            fig = plt.figure(figsize=figsize)
            
            x = uc.get_in_units(x, length_unit)
            disregistry = uc.get_in_units(disregistry, length_unit)
            plt.plot(x, disregistry[:, 0], label='edge disregistry')
            plt.plot(x, disregistry[:, 1], label='normal disregistry')
            plt.plot(x, disregistry[:, 2], label='screw disregistry')
            plt.legend(fontsize='xx-large')
            plt.xlabel(f'x (${length_unit}$)', size='xx-large')
            plt.ylabel(f'disregistry (${length_unit}$', size='xx-large')
            return fig

        else:
            print('x and disregistry must be set/given to plot')

    def E_gsf_surface_plot(self, 
                           x: Optional[npt.ArrayLike] = None,
                           disregistry: Optional[npt.ArrayLike] = None,
                           fmt: str = 'ro-',
                           normalize: Optional[bool] = False,
                           smooth: bool = True, 
                           a1vect: Optional[npt.ArrayLike] = None,
                           a2vect: Optional[npt.ArrayLike] = None,
                           xvect: Optional[npt.ArrayLike] = None,
                           length_unit: str = 'Å',
                           energyperarea_unit: str = 'eV/Å^2',
                           numx: int = 100,
                           numy: int = 100,
                           figsize: Optional[tuple] = None,
                           **kwargs) -> plt.figure:
        """
        Extends the GammaSurface.E_gsf_surface_plot() method to plot the
        disregistry path on top of it.
        
        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates. If x
            or disregistry are not set/given, then the disregistry path will
            not be added.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.  If x
            or disregistry are not set/given, then the disregistry path will
            not be added.
        fmt : str, optional
            The matplotlib.pyplot.plot fmt parameter for the disregistry path
            line, i.e. color, marker and line style options.  Default value is
            'ro-': red with circle markers and solid line.
        normalize : bool, optional
            Flag indicating if axes are Cartesian (False, default) or
            normalized by a1, a2 vectors (True).
        smooth : bool, optional
            If True (default), then plot shows smooth interpolated values.
            If False, plot shows nearest raw data values.
        a1vect : array-like object, optional
            Crystal vector for the a1 vector to use for plotting.  Default
            value of None uses the saved a1vect.
        a2vect : array-like object, optional
            Crystal vector for the a2 vector to use for plotting.  Default
            value of None uses the saved a2vect.
        xvect : array-like object, optional
            Crystal vector to align with the plotting x-axis for 
            non-normalized plots.  If not given, this is taken as the Cartesian
            of a1vect.
        length_unit : str, optional
            The unit of length to display non-normalized axes values in.
            Default value is 'Å'.
        energyperarea_unit : str, optional
            The unit of energy per area to display the stacking fault energies
            in. Default value is 'eV/Å^2'.
        numx : int, optional
            The number of plotting points to use along the x-axis.  Default
            value is 100.
        numy : int, optional
            The number of plotting points to use along the y-axis.  Default
            value is 100.       
        figsize : tuple or None, optional
            The figure's x,y dimensions.  If None (default), the values are
            scaled such that the x,y spacings are approximately equal, and the
            larger of the two values is set to 10.
        **kwargs : dict, optional
            Additional keywords are passed into the underlying 
            matplotlib.pyplot.pcolormesh(). This allows control of such things
            like the colormap (cmap).
            
        Returns
        -------
        matplotlib.pyplot.figure
            The generated figure allowing users to perform additional
            modifications.
        """

        # Generate the surface plot
        fig = self.gamma.E_gsf_surface_plot(normalize=normalize, smooth=smooth, 
                                            a1vect=a1vect, a2vect=a2vect, xvect=xvect,
                                            length_unit=length_unit, energyperarea_unit=energyperarea_unit,
                                            numx=numx, numy=numy, figsize=figsize, **kwargs)
        
        hasvals = True
        # Default values are class properties
        if x is None:
            try:
                x = self.x
            except:
                hasvals = False
        if disregistry is None:
            try:
                disregistry = self.disregistry
            except:
                hasvals = False
        if hasvals:

            # Get xvect direction
            if xvect is None:
                if a1vect is None:
                    a1vect = self.gamma.a1vect
                xvect = np.dot(a1vect, self.gamma.box.vects)

            # Transform disregistry to gamma surface pos
            pos = disregistry.dot(self.transform)
            
            # Transform to x, y plotting coordinates and plot
            x, y = self.gamma.pos_to_xy(pos, xvect=xvect)
            x = uc.get_in_units(x, length_unit)
            y = uc.get_in_units(y, length_unit)
            plt.plot(x, y, fmt)

        return fig

    def E_gsf_vs_x_plot(self, 
                        x: Optional[npt.ArrayLike] = None,
                        disregistry: Optional[npt.ArrayLike] = None,
                        figsize: Optional[tuple] = None, 
                        length_unit: str = 'Å',
                        energyperarea_unit: str = 'eV/Å^2') -> plt.figure:
        """
        Generates a plot of the stacking fault energy, i.e. misfit energy,
        associated with the disregistry values for each x coordinate.
        
        Parameters
        ----------
        x : array-like object, optional
            x-coordinates.  Default value is the stored x-coordinates. If x
            or disregistry are not set/given, then the disregistry path will
            not be added.
        disregistry : array-like object, optional
            (N, 3) shaped array of disregistry vectors at each x-coordinate.
            Default value is the stored disregistry values.  If x
            or disregistry are not set/given, then the disregistry path will
            not be added.
        figsize : tuple or None, optional
            The figure's x,y dimensions.  If None (default), then a figure
            size of (10, 6) will be generated.
        length_unit : str, optional
            The unit of length to display x coordinates in.  Default value is
            'Å'.
        energyperarea_unit : str, optional
            The unit of energy per area to display the stacking fault energies
            in. Default value is 'eV/Å^2'.
            
        Returns
        -------
        matplotlib.pyplot.figure
            The generated figure allowing users to perform additional
            modifications.
        """
        hasvals = True
        # Default values are class properties
        if x is None:
            try:
                x = self.x
            except:
                hasvals = False
        if disregistry is None:
            try:
                disregistry = self.disregistry
            except:
                hasvals = False
        if hasvals:

            if figsize is None:
                figsize=(10,6)

            gsf = self.gamma.E_gsf(pos=disregistry.dot(self.transform))
            
            fig = plt.figure(figsize=figsize)
            plt.plot(uc.get_in_units(x, length_unit),
                     uc.get_in_units(gsf, energyperarea_unit))
            plt.xlabel(f'x-coordinate (${length_unit}$)', size='xx-large')
            plt.ylabel(f'Stacking fault energy (${energyperarea_unit}$)', size='xx-large')
            return fig
        
        else:
            print('x and disregistry must be set/given to check stacking fault energies')


    def load(self, 
             model: Union[str, io.IOBase, DM],
             gamma: Optional[GammaSurface] = None):
        """
        Load solution from a data model.
        
        Parameters
        ----------
        model : str, file-like object or DataModelDict
            The semi-discrete-Peierls-Nabarro data model to load.
        gamma : atomman.defect.GammaSurface, optional
            The gamma surface to use.  If not given, will check to see if the
            content is inside model.
        
        Raises
        ------
        ValueError
            If the gamma surface information is not given and it is not found
            in model.
        """
        
        # Identify model element
        try:
            sdpn = DM(model).find('semidiscrete-variational-Peierls-Nabarro')
        except:
            sdpn = DM(model).find('semi-discrete-Peierls-Nabarro')
        
        # Load calculation parameters
        params = sdpn['parameter']
        try:
            self.__transform = uc.value_unit(params['transform'])
        except: 
            self.__transform = uc.value_unit(params['axes'])
        self.__K_tensor = uc.value_unit(params['K_tensor'])
        self.__burgers = uc.value_unit(params['burgers'])
        self.tau = uc.value_unit(params['tau'])
        self.alpha = uc.value_unit(params['alpha'])
        self.beta = uc.value_unit(params['beta'])
        self.cutofflongrange = uc.value_unit(params['cutofflongrange'])
        self.fullstress = params['fullstress']
        self.cdiffelastic = params['cdiffelastic']
        self.cdiffsurface = params['cdiffsurface']
        self.cdiffstress = params['cdiffstress']
        self.min_method = params['min_method']
        self.min_options = params['min_options']
        self.min_kwargs = params.get('min_kwargs', None)
        
        # Load gamma
        if gamma is None:
            try:
                gamma = GammaSurface(model)
            except:
                raise ValueError('No gamma surface in model or given!')
        elif not isinstance(gamma, GammaSurface):
            gamma = GammaSurface(gamma)
        self.__gamma = gamma
        
        # Load calculation solution
        solution = sdpn['solution']
        self.x = uc.value_unit(solution['x'])
        self.disregistry = uc.value_unit(solution['disregistry'])
    
    def model(self, 
              length_unit: str = 'Å',
              energyperarea_unit: str = 'eV/Å^2',
              pressure_unit: str = 'GPa',
              include_gamma: bool = False) -> DM:
        """
        Generate a data model for the object.
        
        Parameters
        ----------
        length_unit : str, optional
            The unit of length to save values as.  Default is 'Å'.
        energyperarea_unit : str, optional
            The unit of energy per area to save fault energy values as.  Only
            used if the gamma surface information is included.  Default value
            is 'eV/Å^2'.
        pressure_unit : str, optional
            The unit of pressure to save values as.  Default is 'GPa'.
        include_gamma : bool, optional
            Flag indicating if the gamma surface data is to be included.
            Default value is False.
            
        Returns
        -------
        DataModelDict
            The data model containing all input parameters and the current
            disregistry vectors.
        """
        model = DM()
        model['semidiscrete-variational-Peierls-Nabarro'] = sdpn = DM()
        
        sdpn['parameter'] = params = DM()
        params['transform'] = uc.model(self.transform, None)
        params['K_tensor'] = uc.model(self.K_tensor, pressure_unit)
        params['tau'] = uc.model(self.tau, pressure_unit)
        params['alpha'] = uc.model(self.alpha, pressure_unit+'/'+length_unit)
        
        params['beta'] = uc.model(self.beta, pressure_unit+'*'+length_unit)
        params['cdiffelastic'] = self.cdiffelastic
        params['cdiffsurface'] = self.cdiffsurface
        params['cdiffstress'] = self.cdiffstress
        params['cutofflongrange'] = uc.model(self.cutofflongrange, length_unit)
        params['burgers'] = uc.model(self.burgers, length_unit)
        params['fullstress'] = self.fullstress
        params['min_method'] = self.min_method
        params['min_options'] = self.min_options
        if len(self.min_kwargs) > 0:
            params['min_kwargs'] = self.min_kwargs
        
        if include_gamma is True:
            sdpn['generalized-stacking-fault'] = self.gamma.model(
                                                  length_unit=length_unit,
                                                  energyperarea_unit=energyperarea_unit)
        
        sdpn['solution'] = solution = DM()
        solution['x'] = uc.model(self.x, length_unit)
        solution['disregistry'] = uc.model(self.disregistry, length_unit)
        
        return model