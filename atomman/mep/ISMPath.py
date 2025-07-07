# coding: utf-8

# Standard Python imports
from __future__ import annotations
import time
from typing import Optional, Union

from tqdm import trange

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from scipy.interpolate import CubicSpline

from . import BasePath
from ..tools import aslist

class ISMPath(BasePath):
    """
    Class representing an energy path for use with the improved string method.
    """
    
    @property
    def default_timestep(self) -> float:
        """float : The default relaxation timestep"""
        
        # 0.05 * min(0.2, N^-1)
        return 0.05 * np.min([0.2, len(self.coord)**-1])
    
    @property
    def default_tolerance(self) -> float:
        """float : The default relaxation tolerance"""
        
        # max(N^-4, 1e-10)
        return np.max([len(self.coord)**-4, 1e-10])
    
    @property
    def unittangent(self) -> np.ndarray:
        """numpy.NDArray : The tangent vectors along the path at each point."""
        
        τ = np.empty_like(self.coord)
        
        # Calculate all unit difference vectors
        diff = self.coord[1:] - self.coord[:-1]
        diff = (diff.T / np.linalg.norm(diff, axis=-1)).T
        
        # Set end points using only forward, backward
        τ[0] = diff[0]
        τ[-1] = diff[-1]
        
        # Set mid points as average of forward, backward
        τ[1:-1] = diff[:-1] + diff[1:]
        
        # Normalize all to unit vectors
        τ = (τ.T / np.linalg.norm(τ, axis=-1)).T
        
        return τ
    
    def interpolate_path(self,
                         arccoord: npt.ArrayLike) -> ISMPath:
        """
        Uses cubic spline interpolation to interpolate a new path from intermediate
        arc length coordinates along the current path.
        
        Parameters
        ----------
        arccoord: array-like object
            Arc length coordinates along the current path where the new points are
            to be placed.  Values must be in the range [0, self.arccoord].
            
        Returns
        -------
        ISMPath
            A new path with the interpolated coordinates.
        """
        α = self.arccoord
        
        if np.any(arccoord < 0) | np.any(arccoord > α[-1]):
            raise ValueError(f'arccoord values must be in range [0.0, {α[-1]}]')
            
        newcoord = CubicSpline(α, self.coord)(arccoord)
        
        return ISMPath(newcoord, self.energyfxn, gradientfxn=self.gradientfxn,
                    gradientkwargs=self.gradientkwargs)

    def step(self,
             timestep: Optional[float] = None,
             climbindex: Union[int, list, None] = None) -> ISMPath:
        """
        Performs a single string relaxation step.
        
        Parameters
        ----------
        timestep : float, optional
            The size of the timestep to use.  Will use the path's
            default timestep if not given.
        climbindex : int or list, optional
            Indicates the indices of the path points to apply the climb
            algorithm to.  If None, no climbing will be performed.
        
        Returns
        -------
        newpath : ISMPath
            A Path with coordinates evolved forward by one timestep.
        """
        
        # Set default timestep
        if timestep is None:
            timestep = self.default_timestep
        
        # Define the relaxation rate functions
        def rate(coord):
            """Regular motions are along -grad(E)"""
            grad_energy = self.grad_energy(coord)
            return - grad_energy
            
        def climbrate(coord, τ):
            """Climb motions are along -grad(E) + 2 (grad(E) dot τ) τ"""
            grad_energy = self.grad_energy(coord)
            r = - grad_energy + 2 * np.einsum('ij,ij,il->il', grad_energy, τ, τ)
            return r
        
        # Integrate path coordinates using rate
        icoord = self.integratorfxn(rate, self.coord, timestep)
        
        # Integrate the climbing coordinates using climbrate
        if climbindex is not None:
            climbindex = aslist(climbindex)
            
            # Get tangent vectors for initial path
            τ = self.unittangent
            
            # Replace intpath coords for climbindices
            icoord[climbindex] = self.integratorfxn(climbrate,
                                                    self.coord[climbindex],
                                                    timestep,
                                                    τ=τ[climbindex])
        else:
            climbindex = []
        
        # Create intpath from integrated coords
        intpath = ISMPath(icoord, self.energyfxn, self.gradientfxn,
                          self.gradientkwargs)
        
        # Divide full path into segments based on climbindices
        startindices = [0] + aslist(climbindex)
        endindices = aslist(np.asarray(climbindex)+1) + [None]
        
        # Redistribute arclength coords evenly within each segment
        α = intpath.arccoord

        newα = np.empty_like(α)
        for s, e in zip(startindices, endindices):
            
            # Adjust coords in segment to be equally spaced
            subα = α[s:e]
            newα[s:e] = np.linspace(subα[0], subα[-1], len(subα))
        
        # Interpolate new path from the adjusted arclength coords
        newpath = intpath.interpolate_path(newα)
        
        return newpath
    
    def relax(self,
              relaxsteps: int = 0,
              climbsteps: int = 0,
              timestep: Optional[float] = None,
              tolerance: Optional[float] = None,
              climbpoints: int = 1,
              verbose: bool = True) -> ISMPath:
        """
        Perform multiple relaxation and/or climb steps until either the
        maximum coordinate displacement per step drops below a tolerance or
        the maximum number of steps is reached.
        
        Parameters
        ----------
        relaxsteps : int, optional
            The maximum number of relaxation steps to perform.  Default value
            is 0: no relaxation steps.
        climbsteps : int, optional
            The maximum number of climbing steps to perform.  Default value
            is 0: no climbing steps.
        timestep : float, optional
            The size of the timestep to use.  Will use default_timestep if not
            given.
        tolerance : float, optional
            The coordinate displacement tolerance to use.  Will use
            default_tolerance if not given.
        climbpoints : int, optional
            Indicates the maximum number of points to subject the climbing to.
            Default value is 1: i.e. only one maximum is refined.
        verbose : bool, optional
            If True (default), informative statements about the relaxation are
            printed.
        
        Returns
        -------
        newpath : ISMPath
            A Path with coordinates evolved forward from the relaxation.
        """
        
        # Set default timestep
        if timestep is None:
            timestep = self.default_timestep
        
        # Set default tolerance
        if tolerance is None:
            tolerance = self.default_tolerance

        if verbose:
            print(f'timestep =  {timestep}')
            print(f'tolerance = {tolerance}')
            print(flush=True)
        
        # Set current path
        currentpath = self

        # ----------------- Relaxation steps ---------------- #

        if verbose and relaxsteps > 0:
            print('Starting relaxation steps', flush=True)
            itersteps = trange(relaxsteps)
        else:
            itersteps = range(relaxsteps)

        # Perform the relaxation steps
        s = time.time()
        for i in itersteps:
            
            # Step forward
            newpath = currentpath.step(timestep=timestep)
            
            # Calculate max displacement from step
            d = np.linalg.norm(newpath.coord - currentpath.coord, axis=-1).max() / timestep
            
            if verbose:
                itersteps.set_postfix(tolerance=f'{d:.14e}')

            # Update currentpath to newpath
            currentpath = newpath

            # Stop if tolerance reached
            if d < tolerance:
                break
        
        e = time.time()
        #if verbose and relaxsteps > 0:
        #    print(f'Number of relax steps performed: {i+1}')
        #    print(f'Final max displacement: {d}')
        #    print(f'Run time: {e-s} seconds')
        #    print(flush=True)
        
        # --------- Identify climb coordiates(s) -------- #
        
        # Identify the points that are local energy maxima
        energy = currentpath.energy()
        maxmap = np.hstack([False, (energy[1:-1] > energy[:-2]) & (energy[1:-1] > energy[2:]), False])
        climbindex = np.arange(len(maxmap))[maxmap]
        
        if maxmap.sum() > climbpoints:
            climbindex = climbindex[:climbpoints]
        
        if verbose and climbsteps > 0:
            print(f'Starting climbing steps with {len(climbindex)} climbing points', flush=True)
            itersteps = trange(climbsteps)
        else:
            itersteps = range(climbsteps)
        
        # ----------------- Climbing steps ---------------- #
        
        # Perform the climb steps
        s = time.time()
        for i in itersteps:
            
            # Step forward
            newpath = currentpath.step(timestep=timestep, climbindex=climbindex)
            
            # Calculate max displacement
            d = np.linalg.norm(newpath.coord - currentpath.coord, axis=-1).max() / timestep
            
            if verbose:
                itersteps.set_postfix(tolerance=f'{d:.14e}')
            
            # Update currentpath to newpath
            currentpath = newpath
            
            # Stop if tolerance reached
            if d < tolerance:
                break
        
        e = time.time()
        if verbose and climbsteps > 0:
            #print(f'Number of climb steps performed: {i+1}')
            #print(f'Final max displacement: {d}')
            print(f'Run time: {e-s} seconds')
            print(f'Max energy before climb = {energy.max()}', flush=True)
            print(f'Max energy after climb = {currentpath.energy().max()}')
            print(flush=True)
        
        return currentpath