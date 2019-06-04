# Standard Python imports
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import os
from copy import deepcopy

# http://www.numpy.org/
import numpy as np

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

from ..tools import axes_check, vect_angle, miller

class StackingFault(object):
    def __init__(self, system, a1vect=None, a2vect=None,
                 ucellbox=None, transform=None, cutboxvector='c',
                 faultposrel=None, faultposcart=None, vacuumwidth=None):
        """
        Initializes an object for generating generalized stacking fault
        atomic configurations.

        Parameters
        ----------
        system : atomman.System
            The atomic system to generate stacking fault configurations from.
            The system should be oriented such that the stacking fault plane
            is parallel to either the xy-plane, the xz-plane or the yz-plane.
        a1vect : array-like object, optional
            A slip vector within the slip plane.  Depending on if ucellbox and
            transform are given, this can be either a Miller crystal vector or
            a Cartesian vector relative to the supplied system.  If a1vect is
            not given and a2vect is, then a1vect is set to [0,0,0].
        a2vect : array-like object, optional
            A slip vector within the slip plane.  Depending on if ucellbox and
            transform are given, this can be either a Miller crystal vector or
            a Cartesian vector relative to the supplied system.  If a2vect is
            not given and a1vect is, then a2vect is set to [0,0,0].
        ucellbox : atomman.Box, optional
            If ucellbox is given, then the a1vect and a2vect slip vectors are
            taken as Miller crystal vectors relative to ucellbox.  If not
            given, then the slip vectors are taken as Cartesian vectors.
        transform : array-like object, optional
            A transformation tensor to apply to the a1vect and a2vect slip
            vectors.  This is needed if system is oriented differently than
            ucellbox, i.e. system is rotated.
        cutboxvector : str, optional
            Defines which of the three box vectors gets cut so that the slip
            plane can be added.  'a' indicates the slip plane is parallel to
            the yz-plane, 'b' the xz-plane and 'c' the xy-plane.  Note that 
            the specified cutboxvector is the only box vector allowed to have a
            component out of the slip plane.
        faultposrel : float, optional
            The position to place the slip plane within the system given as a
            relative coordinate along the out-of-plane direction.  faultposrel
            and faultposcart cannot both be given.  Default value is 0.5 if 
            faultposcart is also not given.
        faultposcart : float, optional
            The position to place the slip plane within the system given as a
            Cartesian coordinate along the out-of-plane direction.  faultposrel
            and faultposcart cannot both be given.
        vacuumwidth : float, optional
            If given, the out-of-plane dimension of system is increased by this
            amount to create a vacuum region.  Useful if the generated
            configurations are to be evaluated with DFT.
        """
        
        self.set(system, a1vect=a1vect, a2vect=a2vect, ucellbox=ucellbox,
                    transform=transform, cutboxvector=cutboxvector,
                    faultposrel=faultposrel, faultposcart=faultposcart,
                    vacuumwidth=vacuumwidth)
        

    def set(self, system, a1vect=None, a2vect=None, ucellbox=None, transform=None,
            cutboxvector='c', faultposrel=None, faultposcart=None,
            vacuumwidth=None):
        """
        Sets generation parameters associated with a stacking fault.

        Parameters
        ----------
        system : atomman.System
            The atomic system to generate stacking fault configurations from.
            The system should be oriented such that the stacking fault plane
            is parallel to either the xy-plane, the xz-plane or the yz-plane.
        a1vect : array-like object, optional
            A slip vector within the slip plane.  Depending on if ucellbox and
            transform are given, this can be either a Miller crystal vector or
            a Cartesian vector relative to the supplied system.  If a1vect is
            not given and a2vect is, then a1vect is set to [0,0,0].
        a2vect : array-like object, optional
            A slip vector within the slip plane.  Depending on if ucellbox and
            transform are given, this can be either a Miller crystal vector or
            a Cartesian vector relative to the supplied system.  If a2vect is
            not given and a1vect is, then a2vect is set to [0,0,0].
        ucellbox : atomman.Box, optional
            If ucellbox is given, then the a1vect and a2vect slip vectors are
            taken as Miller crystal vectors relative to ucellbox.  If ucellbox
            is a standard hexagonal cell, the slip vectors can alternatively be
            given as Miller-Bravais crystal vectors. If not given, then the
            slip vectors are taken as Cartesian vectors.
        transform : array-like object, optional
            A transformation tensor to apply to the a1vect and a2vect slip
            vectors.  This is needed if system is oriented differently than
            ucellbox, i.e. system is rotated.
        cutboxvector : str, optional
            Defines which of the three box vectors gets cut so that the slip
            plane can be added.  'a' indicates the slip plane is parallel to
            the yz-plane, 'b' the xz-plane and 'c' the xy-plane.  Note that 
            the specified cutboxvector is the only box vector allowed to have a
            component out of the slip plane.
        faultposrel : float, optional
            The position to place the slip plane within the system given as a
            relative coordinate along the out-of-plane direction.  faultposrel
            and faultposcart cannot both be given.  Default value is 0.5 if 
            faultposcart is also not given.
        faultposcart : float, optional
            The position to place the slip plane within the system given as a
            Cartesian coordinate along the out-of-plane direction.  faultposrel
            and faultposcart cannot both be given.
        """
        # Copy system for safe manipulation
        system = deepcopy(system)
        system.wrap()

        # Check fault pos values
        if faultposcart is None and faultposrel is None:
            faultposrel = 0.5 
        elif faultposcart is not None and faultposrel is not None:
            raise ValueError('faultposrel and faultposcart cannot both be given')

        # Set cutindex and faultarea based on cutboxvector
        if cutboxvector == 'a':
            if system.box.bvect[0] != 0.0 or system.box.cvect[0] != 0.0:
                raise ValueError("box bvect and cvect cannot have x component for cutboxvector='a'")
            cutindex = 0
            faultarea = np.linalg.norm(np.cross(system.box.bvect, system.box.cvect))
        
        elif cutboxvector == 'b':
            if system.box.avect[1] != 0.0 or system.box.cvect[1] != 0.0:
                raise ValueError("box avect and cvect cannot have y component for cutboxvector='b'")
            cutindex = 1
            faultarea = np.linalg.norm(np.cross(system.box.avect, system.box.cvect))

        elif cutboxvector == 'c':
            if system.box.avect[2] != 0.0 or system.box.bvect[2] != 0.0:
                raise ValueError("box avect and bvect cannot have z component for cutboxvector='c'")
            cutindex = 2
            faultarea = np.linalg.norm(np.cross(system.box.avect, system.box.bvect))

        else: 
            raise ValueError('Invalid cutboxvector')

        # Define out of plane unit vector 
        ovect = np.zeros(3)
        ovect[cutindex] = 1.0

        # Set system's pbc
        system.pbc = [True, True, True]
        system.pbc[cutindex] = False

        # Add vacuum to non-periodic direction
        if vacuumwidth is not None:
            if vacuumwidth < 0:
                raise ValueError('vacuumwidth must be positive')
            newvects = system.box.vects
            newvects[cutindex, cutindex] += vacuumwidth
            neworigin = system.box.origin - ovect * vacuumwidth / 2 
            system.box_set(vects=newvects, origin=neworigin)

        # Identify atoms above fault plane position
        if faultposcart is None:
            faultposcart = system.box.origin[cutindex] + system.box.vects[cutindex, cutindex] * faultposrel            
        abovefault = system.atoms.pos[:, cutindex] > (faultposcart)

        # Set default a1vect or a2vect if needed
        if a1vect is None and a2vect is None:
            raise ValueError('At least one of a1vect and a2vect must be given')
        elif a1vect is None:
            a1vect = np.zeros(3)
        elif a2vect is None:
            a2vect = np.zeros(3)

        # Convert crystal vects to Cartesian
        if ucellbox is not None:
            a1vectcart = miller.vector_crystal_to_cartesian(a1vect, ucellbox)
            a2vectcart = miller.vector_crystal_to_cartesian(a2vect, ucellbox)
        else:
            a1vectcart = a1vect
            a2vectcart = a2vect
        
        # Transform vects to match system
        if transform is not None:
            transform = axes_check(transform)
            a1vectcart = transform.dot(a1vectcart)
            a2vectcart = transform.dot(a2vectcart)

        # Check shift vectors values
        if not np.isclose(a1vectcart[cutindex], 0.0):
            raise ValueError('a1vect not in slip plane')
        if not np.isclose(a2vectcart[cutindex], 0.0):
            raise ValueError('a2vect not in slip plane')
        if not np.allclose(a1vectcart, np.zeros(3)) and not np.allclose(a2vectcart, np.zeros(3)):
            angle = vect_angle(a1vectcart, a2vectcart)
            if np.isclose(angle, 0.0) or np.isclose(angle, 180.0):
                raise ValueError('a1vect and a2vect cannot be parallel')
    
        # Save properties
        self.__system = system
        self.__a1vect = a1vect
        self.__a2vect = a2vect
        self.__a1vectcart = a1vectcart
        self.__a2vectcart = a2vectcart
        self.__cutboxvector = cutboxvector
        self.__cutindex = cutindex
        self.__faultposcart = faultposcart
        self.__abovefault = abovefault
        self.__faultarea = faultarea
        self.__ucellbox = ucellbox
        self.__transform = transform

    @property
    def system(self):
        """atomman.System : The unshifted configuration."""
        return self.__system

    @property
    def a1vect(self):
        """numpy.ndarray : One of the two shift vectors as given."""
        return self.__a1vect

    @property
    def a2vect(self):
        """numpy.ndarray : One of the two shift vectors as given."""
        return self.__a2vect

    @property
    def a1vectcart(self):
        """numpy.ndarray : One of the two shift vectors in Cartesian relative to system."""
        return self.__a1vectcart
    
    @property
    def a2vectcart(self):
        """numpy.ndarray : One of the two shift vectors in Cartesian relative to system."""
        return self.__a2vectcart

    @property
    def cutboxvector(self):
        """str : The box vector that is cut to add vacuum and slip plane."""
        return self.__cutboxvector

    @property
    def cutindex(self):
        """int : The Cartesian index for the cut direction."""
        return self.__cutindex

    @property
    def faultposcart(self):
        """float : The Cartesian position of the slip plane."""
        return self.__faultposcart

    @property
    def abovefault(self):
        """list : Indices of all atoms in System above the slip plane."""
        return self.__abovefault
    
    @property
    def faultarea(self):
        """float : The area of the slip plane."""
        return self.__faultarea

    @property
    def ucellbox(self):
        """atomman.Box or None : The supplied reference cell box."""
        return self.__ucellbox

    @property
    def transform(self):
        """numpy.ndarray or None : The supplied transformation matrix."""
        return self.__transform
    
    def fault(self, a1=None, a2=None, outofplane=None, faultshift=None, minimum_r=None):
        """
        Generates a fault configuration by displacing all atoms above the slip
        plane.

        Parameters
        ----------
        a1 : float, optional
            The fractional coordinate of a1vect to shift by. 
            Default value is 0.0.
        a2 : float, optional
            The fractional coordinate of a2vect to shift by. 
            Default value is 0.0. 
        outofplane : float, optional
            An out-of-plane shift, given in absolute units.
            Default value is 0.0.
        faultshift : array-like object, optional
            The full shifting vector to displace the atoms above the slip
            plane by.  Cannot be given with a1, a2, or outofplane.
        minimum_r : float, optional
            Specifies the minimum allowed interatomic spacing across the slip
            plane.  If any sets of atoms are closer than this value then the
            outofplane shift is increased.  Default value is None, which
            performs no adjustment.

        Returns
        -------
        atomman.System
            The atomic configuration with stacking fault shift
        """
        # Define out of plane unit vector 
        ovect = np.zeros(3)
        ovect[self.cutindex] = 1.0
        
        # Identify the two non-cut indices
        inindex = []
        for i in range(3):
            if i != self.cutindex:
                inindex.append(i)
        
        # Calculate faultshift
        if a1 is not None or a2 is not None or outofplane is not None:
            if faultshift is not None:
                raise ValueError('a1, a2, outofplane cannot be given with faultshift')
            if a1 is None:
                a1 = 0.0
            if a2 is None:
                a2 = 0.0
            if outofplane is None:
                outofplane = 0.0
            faultshift = a1 * self.a1vectcart + a2 * self.a2vectcart + outofplane * ovect
        
        # Set default faultshift
        elif faultshift is None:
            faultshift = np.array([0.0, 0.0, 0.0])

        # Shift atoms above the fault by faultshift
        sfsystem = deepcopy(self.system)
        sfsystem.atoms.pos[self.abovefault] += faultshift
        sfsystem.wrap()
        
        # Add additional outofplane shift if necessary
        if minimum_r is not None:
            # Get all atoms within minimum_r of fault position
            top_pos = sfsystem.atoms.pos[self.abovefault]
            top_pos = top_pos[top_pos[:, self.cutindex] <= self.faultposcart + minimum_r] 
            bot_pos = sfsystem.atoms.pos[~self.abovefault]
            bot_pos = bot_pos[bot_pos[:, self.cutindex] >= self.faultposcart - minimum_r]
            if top_pos.shape[0] > 0 and bot_pos.shape[0] > 0:
                dmag_min = minimum_r
                dvect_min = None
                
                for i in range(top_pos.shape[0]):
                    dvect = sfsystem.dvect(bot_pos, top_pos[i])
                    if dvect.shape == (3,):
                        dvect = dvect.reshape(1,3)
                    dmag = np.linalg.norm(dvect, axis=1)
                    i = np.argmin(dmag)
                    if dmag[i] < dmag_min:
                        dmag_min = dmag[i]
                        dvect_min = dvect[i]
                
                if dvect_min is not None:
                    
                    new = (minimum_r**2 - dvect_min[inindex[0]]**2 - dvect_min[inindex[1]]**2)**0.5
                    outofplane = new - dvect_min[self.cutindex]
                    faultshift = outofplane * ovect
                    sfsystem.atoms.pos[self.abovefault] += faultshift
                    sfsystem.wrap()
        
        return sfsystem
    
    def iterfaultmap(self, num_a1=None, num_a2=None, outofplane=None, minimum_r=None):
        """
        Iterates over generalized stacking fault configurations associated
        with a 2D map of equally spaced a1, a2 coordinates.

        Parameters
        ----------
        num_a1 : int
            The number of a1 values to generate systems for.  
            Default value is 1 (only generate for a1=0.0).
        num_a2 : int
            The number of a2 values to generate systems for.  
            Default value is 1 (only generate for a2=0.0).
        outofplane : float, optional
            An out-of-plane shift, given in absolute units.
            Default value is 0.0.
        minimum_r : float, optional
            Specifies the minimum allowed interatomic spacing across the slip
            plane.  If any sets of atoms are closer than this value then the
            outofplane shift is increased.  Default value is None, which
            performs no adjustment.
        
        Yields
        ------
        a1 : float
            The a1 fractional coordinate of a1vect. 
        a2 : float
            The a2 fractional coordinate of a2vect. 
        atomman.System
            The fault configuration associated with the a1, a2 shift.
        """
        if num_a1 is None:
            num_a1 = 1
        if num_a2 is None:
            num_a2 = 1

        # Construct mesh of regular points
        a1s, a2s = np.meshgrid(np.linspace(0, 1, num_a1, endpoint=False),
                               np.linspace(0, 1, num_a2, endpoint=False))

        for a1, a2 in zip(a1s.flat, a2s.flat):
            yield a1, a2, self.fault(a1=a1, a2=a2, outofplane=outofplane,
                                     minimum_r=minimum_r)