# coding: utf-8
# Standard Python libraries
from __future__ import annotations
import io
from collections import OrderedDict
from copy import deepcopy
import warnings
from typing import Any, Optional, Union, Tuple

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# https://pandas.pydata.org/
import pandas as pd

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
import atomman.unitconvert as uc
from . import Atoms, Box, dvect, dmag, NeighborList
from ..lammps import normalize as lmp_normalize
from ..tools import indexstr, miller, ishexagonal, aslist
from .. import dump

class System(object):
    """
    A representation of an atomic system.  This combines Atoms with Box and
    adds methods and attributes involving both.
    """
    
    class _AtomsIndexer(object):
        """Internal class for setitem / getitem acting on atoms"""
        def __init__(self, host: System):
            self.__host = host
            
        def __getitem__(self,
                        index: Union[int, list, slice]) -> System:
            """Index getting of Atoms that operates on System."""
            host = self.__host
            return System(atoms=host.atoms[index], box=host.box, pbc=host.pbc,
                          symbols=host.symbols)
            
        def __setitem__(self,
                        index: Union[int, list, slice],
                        value: Any):
            """Index setting of Atoms that operates on System."""
            host = self.__host
            if isinstance(value, Atoms):
                host.atoms[index] = value
            elif isinstance(value, System):
                try:
                    assert np.allclose(host.box.vects, value.box.vects)
                    assert np.allclose(host.box.origin, value.box.origin)
                except:
                    warnings.warn('Atom assignment between two Systems with different boxes', UserWarning)
                host.atoms[index] = value.atoms
            else:
                raise ValueError('Can only set using Atoms or System objects')
    
    def __init__(self,
                 atoms: Optional[Atoms] = None,
                 box: Optional[Box] = None,
                 pbc: Optional[Tuple[bool, bool, bool]] = None,
                 scale: bool = False,
                 symbols: Union[str, list, None] = None,
                 masses: Union[float, list, None] = None,
                 model: Union[str, io.IOBase, DM, None] = None,
                 safecopy: bool = False):
        """
        Initialize a System by joining an am.Atoms and am.Box instance.
        
        Parameters
        ----------
        atoms : atomman.Atoms, optional
            The underlying Atoms object to build system around.
        box : atomman.Box, optional
            The underlying box object to build system around.
        pbc : tuple or list of bool, optional
            Indicates which of the dimensions related to the three box vectors
            are periodic.  Default value is (True, True, True).
        scale : bool, optional
            If True, atoms.pos will be scaled relative to the box.  Default
            value is False.
        symbols : tuple, optional
            A list of the element symbols for each atom atype.  If len(symbols)
            is less than natypes, then missing values will be set to None.
            Default sets list with all None values.
        masses : tuple, optional
            A list of the masses for each atom atype.  If len(symbols) is less
            than natypes, then missing values will be set to None.  Default
            sets list with all None values.
        model : str or DataModelDict, optional
            File path or content of a JSON/XML data model containing all
            system information.  Cannot be given with atoms, box or scale.
        safecopy : bool, optional
            Flag indicating if values are to be copied before setting.  For
            values given as objects, direct setting (False, default) may result
            in the System pointing to the original object.  Using safecopy=True
            deep copies the objects before setting to avoid this.  Note that
            safecopy=True may be considerably slower for large numbers of atoms
            and/or properties.
        
        Returns
        =======
        System
            The System object.
        """
        # Check for model
        if model is not None:
            try:
                assert atoms is None
                assert box is None
                assert scale is False
            except:
                raise ValueError('model cannot be given with atoms, box or scale parameters')
            
            # Load data model
            model = DM(model).find('atomic-system')

            # Extract values from model
            box = Box(model=model)
            atoms = Atoms(model=model)
            if pbc is None:
                pbc = model['periodic-boundary-condition']
            if symbols is None:
                symbols = tuple(model.aslist('atom-type-symbol'))
            if masses is None:
                masses = tuple(model.aslist('atom-type-mass'))

        # Interpret given/missing parameters
        else:
            # Set default atoms or deepcopy
            if atoms is None:
                atoms = Atoms()
            elif safecopy:
                atoms = deepcopy(atoms)

            # Set default box or deepcopy
            if box is None:
                box = Box()
            elif safecopy:
                box = deepcopy(box)
            
            # Set default pbc
            if pbc is None:
                pbc = (True, True, True)
            
            # Set default masses
            if masses is None:
                masses = []
            else:
                masses = aslist(masses)
            
            # Set default symbols
            if symbols is None:
                symbols = [None for i in range(len(masses))]
            else:
                symbols = aslist(symbols)
            
            # Check data types
            if not isinstance(atoms, Atoms):
                raise TypeError('Invalid atoms type')
            if not isinstance(box, Box):
                raise TypeError('Invalid box type')
            if not isinstance(scale, bool):
                raise TypeError('Invalid scale type')
        
        # Set properties
        self.__atoms = atoms
        self.__box = box
        self.pbc = pbc
        self.__transformation = np.identity(3)

        self.symbols = symbols
        self.masses = masses
        
        # Scale pos if needed
        if scale is True:
            self.atoms_prop('pos', value=self.atoms.pos, scale=True)
        
        # Scale model properties if needed
        if model is not None:
            for prop in model['atoms'].aslist('property'):
                if prop['data'].get('unit', None) == 'scaled':
                    self.atoms.view[prop['name']] = self.box.position_relative_to_cartesian(self.atoms.view[prop['name']]) 
        
        # Set atoms indexer
        self.__atoms_ix = System._AtomsIndexer(self)
    
    def __str__(self) -> str:
        """str : The string representation of a system."""
        return '\n'.join([str(self.box),
                          'natoms = ' + str(self.natoms),
                          'natypes = ' + str(self.natypes),
                          'symbols = ' + str(self.symbols),
                          'pbc = ' + str(self.pbc),
                          str(self.atoms)])
    
    def __len__(self) -> int:
        return self.natoms
    
    @property
    def atoms(self) -> Atoms:
        """atomman.Atoms : underlying Atoms object."""
        return self.__atoms
    
    @property
    def natoms(self) -> int:
        """int : The number of atoms in the Atoms class."""
        return self.__atoms.natoms
    
    @property
    def atypes(self) -> tuple:
        """tuple : List of int atom types."""
        return tuple(range(1, self.natypes+1))
    
    @property
    def natypes(self) -> int:
        """int : The number of atom types/symbols."""
        
        # Check if more symbols than atypes
        try:
            nsymbols = len(self.symbols)
        except:
            nsymbols = 0
        if nsymbols > self.__atoms.natypes:
            return len(self.symbols)
        else:
            return self.__atoms.natypes
    
    @property
    def box(self) -> Box:
        """atommman.Box : underlying Box object."""
        return self.__box
    
    @property
    def pbc(self) -> np.ndarray:
        """numpy.ndarray of bool : The periodic boundary condition settings."""
        return self.__pbc
        
    @pbc.setter
    def pbc(self, value: npt.ArrayLike):
        pbc = np.asarray(value, dtype=bool)
        assert pbc.shape == (3,), 'invalid pbc entry' 
        self.__pbc = pbc
    
    @property
    def atoms_ix(self) -> _AtomsIndexer:
        """Indexer for index slicing of Systems by per-atom properties."""
        return self.__atoms_ix

    @property
    def symbols(self) -> tuple:
        """tuple : The element model symbols associated with each atype."""
        
        # Fill in missing values
        if len(self.__symbols) < self.__atoms.natypes:
            self.symbols = self.__symbols
        return self.__symbols
    
    @symbols.setter
    def symbols(self, value: Union[str, list]):
        
        # Make value list if needed
        value = aslist(value)
        
        # Fill in missing values
        if len(value) < self.__atoms.natypes:
            newvalue = [None for x in range(self.__atoms.natypes)]
            for i in range(len(value)):
                newvalue[i] = value[i]
            value = newvalue

        self.__symbols = tuple(value)

    @property
    def masses(self) -> tuple:
        """tuple : The masses associated with each atype if given."""
        
        # Fill in missing values
        if len(self.__masses) < self.natypes:
            self.masses = self.__masses
        return self.__masses
    
    @masses.setter
    def masses(self, value: Union[float, list]):
        
        # Make value list if needed
        value = aslist(value)
        
        # Convert non-None values to float
        for i in range(len(value)):
            if value[i] is not None:
                value[i] = float(value[i])

        # Fill in missing values
        if len(value) < self.natypes:
            newvalue = [None for x in range(self.natypes)]
            for i in range(len(value)):
                newvalue[i] = value[i]
            value = newvalue
        elif len(value) > self.natypes:
            raise ValueError('More masses than atom types given. Either change atype values or symbols first.')

        self.__masses = tuple(value)
    
    @property
    def composition(self) -> Optional[str]:
        """
        The system's reduced and sorted symbols composition.
        
        Returns
        -------
        str or None
            The system's reduced and sorted symbols composition.
            Will return None if any symbols are missing
        """
        
        # Build total count of each symbol
        sym_dict = {}
        for i in range(self.natypes):
            count = np.sum(self.atoms.atype == i+1)
            if count > 0:
                symbol = self.symbols[i]
                
                if symbol is None:
                    return None
                
                if symbol in sym_dict:
                    sym_dict[symbol] += count
                else:
                    sym_dict[symbol] = count

        # Find greatest common denominator
        gcd = np.gcd.reduce(list(sym_dict.values())) # pylint: disable=no-member
        
        # Sort symbols and reduce counts
        composition =''
        for symbol in sorted(sym_dict):
            count = sym_dict[symbol] // gcd
            if sym_dict[symbol] > 0:
                composition += symbol
                if count != 1:
                    composition += str(count)
        
        return composition
    
    def atoms_prop(self,
                   key: Optional[str] = None,
                   index: Union[int, list, slice, None] = None,
                   value: Optional[Any] = None,
                   a_id: Optional[int] = None,
                   scale: bool = False
                   ) -> Union[list, Atoms, np.ndarray, None]:
        """
        Extends Atoms.prop() by adding a scale argument.
        
        Parameters
        ----------
        key : str, optional
            Per-atom property name.
        index : int, list, slice, optional
            Index of atoms.
        value : any, optional
            Property values to assign.
        a_id : int, optional
            Integer atom index.  Left in for backwards compatibility.
        scale : bool, optional
            Flag indicating if values should be scaled/unscaled. If True:
            - Values being retrieved are scaled from absolute Cartesian to box relative vectors.
            - Values being set are unscaled from box relative to absolute Cartesian vectors.
            Default value is False.
        
        Returns
        -------
        list
            If no parameters given, returns a list of all assigned property
            keys.
        atomman.Atoms
            If index or a_id is given without value or key, returns a new
            Atoms instance for the specified atom indices.
        numpy.ndarray
            If key (and index/a_id) is given without value, returns a copy of
            the data associated with that property key.
        """
        
        # Check that scale is bool
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Call atoms.prop() for scale is False
        if scale is False:
            if value is None:
                return self.atoms.prop(key=key, index=index, a_id=a_id)
            else:
                self.atoms.prop(key=key, index=index, value=value, a_id=a_id)
        
        # Handle property scaling if scale is True
        else:
            # Handle a_id
            if a_id is not None:
                if index is not None:
                    raise ValueError('a_id and index cannot both be given')
                index = a_id
            
            # Get values if value is None
            if value is None:
                # If no key, return copy of atoms (all or a slice) with scaled pos
                if key is None:
                    if index is None:
                        newatoms = deepcopy(self.atoms)
                    else:
                        newatoms = deepcopy(self.atoms[index])
                    newatoms.pos = self.box.position_cartesian_to_relative(newatoms.pos)
                    return newatoms
                
                # If key is given, return scaled property values
                else:
                    if index is None:
                        value = self.atoms.view[key]
                    else:
                        value = self.atoms.view[key][index]
                    return self.box.position_cartesian_to_relative(value)
            
            # Set values if value is given
            else:
                # If no key, unscale pos of Atoms value and set to atoms (or a slice)
                if key is None:
                    if not isinstance(value, Atoms):
                        raise TypeError('If key is None, value must be instance of atomman.Atoms')
                    
                    value.pos = self.box.position_relative_to_cartesian(value.pos)
                    if index is None:
                        self.atoms[:] = value
                    else:
                        self.atoms[index] = value
                
                # If key is given, unscale value and set to property
                else:
                    value = self.box.position_relative_to_cartesian(value)
                    if index is None:
                        self.atoms.view[key] = value
                    else:
                        self.atoms.view[key][index] = value
    
    def atoms_df(self,
                 scale: bool = False) -> pd.DataFrame:
        """
        Extends Atoms.df() by adding a scale argument.
        
        Parameters
        ----------
        scale : bool or list, optional
            Indicates if/which per-atom properties are to be scaled to box
            relative values.  If False (default), no properties will be scaled.
            If True, pos will be scaled.  Multiple properties can be scaled by
            providing a list of the property names to scale.
            
        Returns
        -------
        pandas.DataFrame
            Tabulated DataFrame version of the atomic data.
        """
        # Handle scale values
        if scale is True:
            scale = ['pos']
        elif scale is False:
            scale = []
        elif not isinstance(scale, list):
            scale = [scale]
        
        # Initialize new dictionary of values
        values = OrderedDict()
        for key in self.atoms.view.keys():
            
            value = self.atoms.view[key]
            if key in scale:
                value = self.box.position_cartesian_to_relative(value)
            
            # Flatten multidimensional arrays
            for index, istr in indexstr(self.atoms.view[key].shape[1:]):
                newkey = key + istr
                
                # Copy values over
                if index == ():
                    values[newkey] = value
                else:
                    values[newkey] = value[(Ellipsis, ) + index]
        
        # Return DataFrame
        return pd.DataFrame(values)
    
    def atoms_extend(self,
                     value: Union[Atoms, int],
                     scale: bool = False,
                     symbols: Optional[list] = None,
                     safecopy: bool = False) -> System:
        """
        Extends Atoms.extend() to the System level by adding scale and symbols
        arguments.
        
        Parameters
        ----------
        value : atomman.Atoms or int
            An int value will result in the atoms object being extended by
            that number of atoms, with all per-atom properties having default
            values (atype = 1, everything else = 0).  For an Atoms value, the
            current atoms list will be extended by the correct number of atoms
            and all per-atom properties in value will be copied over.  Any
            properties defined in one Atoms object and not the other will be
            set to default values.
        scale : bool, optional
            Flag indicating if position values in a supplied Atoms value are to
            be taken as absolute Cartesian (False, default) or in scaled box
            relative units (True).
        symbols : tuple, list or None, optional
            Allows for the system's symbols list to be updated.  If not given,
            will use the current object's symbols.
        safecopy : bool, optional
            Flag indicating if values are to be copied before setting.  If 
            False (default), underlying objects may be shared between the new
            system and the current system and input parameters.  If True, atoms
            and box will be deepcopied before setting.
            Note that safecopy=True may be considerably slower for large
            numbers of atoms and/or properties.
        
        Returns
        -------
        atomman.System
            A new System object with Atoms extended to contain all atoms and
            properties of the current object plus the additional atoms.  The
            current System object's box and pbc (and symbols if not specified)
            will be copied over.
        """
        # Copy value if safecopy
        if safecopy:
            value = deepcopy(value)
        
        # scale only makes sense for Atoms values
        if scale is True and not isinstance(value, Atoms):
            raise ValueError('scale can only be True for Atoms values')        
        
        # Handle symbols parameter
        if symbols is None:
            symbols = self.symbols
        
        # Handle box
        if safecopy:
            box = deepcopy(self.box)
        else:
            box = self.box
        
        # Call atoms.extend to generate new atoms
        atoms = self.atoms.extend(value)
        
        # Unscale pos from Atoms value if needed
        if scale:
            atoms.pos[value.natoms:] = self.box.position_relative_to_cartesian(value.pos)

        # Generate and return new System
        return System(atoms=atoms, box=box, pbc=self.pbc, symbols=symbols)
    
    def box_set(self, **kwargs):
        """
        Extends box.set() with a scale argument.
        
        Parameters
        ----------
        scale : bool
            If True, the scaled (box-relative) positions remain unchanged.
            Default value is False (absolute positions unchanged).
        other kwargs : any
            Any other kwargs allowed by box.set().
        """
        
        # Pop scale
        scale = kwargs.pop('scale', False)
        if not isinstance(scale, bool):
            raise TypeError('Invalid scale type')
        
        # Hold scaled positions constant
        if scale is True:
            spos = self.atoms_prop('pos', scale=True)
            self.box.set(**kwargs)
            self.atoms_prop('pos', value=spos, scale=True)
        
        # Call box.set without scaling
        else:
            self.box.set(**kwargs)
    
    def scale(self, value: npt.ArrayLike) -> np.ndarray:
        """
        Scales 3D vectors from absolute Cartesian coordinates to relative box
        coordinates.
        
        Parameters
        ----------
        value : array-like object
            Absolute Cartesian coordinates to scale.

        Returns
        -------
        numpy.ndarray
            The relative box coordinates associated with the given values.
        """
        warnmsg = "System.scale() will likely be depreciated in the next big version update "
        warnmsg += "as the method name is not informative.  It is being replaced by the "
        warnmsg += "equivalent Box.position_cartesian_to_relative() method."
        warnings.warn(warnmsg, PendingDeprecationWarning)
        
        # Retrieve parameters
        #value = np.asarray(value, dtype=float)
        #vects = self.box.vects
        #inverse = np.linalg.inv(vects)
        #origin = self.box.origin
        
        # Convert
        #return (value - origin).dot(inverse)

        return self.box.position_cartesian_to_relative(value)

    def unscale(self, value: npt.ArrayLike) -> np.ndarray:
        """
        Unscales 3D vectors from relative box coordinates to absolute
        Cartesian coordinates.
        
        Parameters
        ----------
        value : numpy.ndarray
            Relative box coordinates to unscale.

        Returns
        -------
        numpy.ndarray
            Absolute Cartesian coordinates associated with the given values.
        """
        warnmsg = "System.unscale() will likely be depreciated in the next big version update "
        warnmsg += "as the method name is not informative.  It is being replaced by the "
        warnmsg += "equivalent Box.position_relative_to_cartesian() method."
        warnings.warn(warnmsg, PendingDeprecationWarning)
        
        # Retrieve parameters
        #value = np.asarray(value, dtype=float)
        #vects = self.box.vects
        #origin = self.box.origin
        
        # Convert
        #return value.dot(vects) + origin

        return self.box.position_relative_to_cartesian(value)
        
    def wrap(self,
             return_imageflags: bool = False) -> Optional[np.ndarray]:
        """
        Wrap atoms around periodic boundaries and extend non-periodic
        boundaries such that all atoms are within the box.

        Parameters
        ----------
        return_imageflags : bool, optional
            If True, an array of which image the atom was originally in is
            returned.

        Returns
        -------
        numpy.NDarray of int
            The imageflags array - only returned if return_imageflags = True.
        """
        
        # mins and maxs are box dimensions relative to box vectors, i.e 0 to 1
        mins = np.array([0.0, 0.0, 0.0])
        maxs = np.array([1.0, 1.0, 1.0])
        
        # Retrieve scaled pos
        spos = self.atoms_prop('pos', scale=True)

        # Initialize wrapflags
        imageflags = np.zeros_like(spos, dtype=int)
        
        # Loop over three pbc directions
        for i in range(3):
            
            # Count wraps across periodic boundaries
            if self.pbc[i]:
                imageflags[:, i] = np.floor(spos[:, i]) # pylint: disable=unsupported-assignment-operation
            
            # Shift min and max to encompass atoms across non-periodic bounds
            else:
                min = spos[:, i].min()
                max = spos[:, i].max()
                if min <= mins[i]: 
                    mins[i] = min - 0.001
                if max >= maxs[i]: 
                    maxs[i] = max + 0.001
            
        # Wrap atoms across periodic boundaries
        spos -= imageflags
        
        # Unscale spos and save to pos
        self.atoms_prop('pos', value=spos, scale=True)
        
        # Modify box vectors and origin by new min and max
        origin = self.box.origin + mins.dot(self.box.vects) 
        avect = self.box.avect * (maxs[0] - mins[0])
        bvect = self.box.bvect * (maxs[1] - mins[1])
        cvect = self.box.cvect * (maxs[2] - mins[2])
        self.box_set(avect=avect, bvect=bvect, cvect=cvect, origin=origin)

        if return_imageflags:
            return imageflags
    
    def dvect(self,
              pos_0: Union[int, list, slice, npt.ArrayLike],
              pos_1: Union[int, list, slice, npt.ArrayLike]
              ) -> np.ndarray:
        """
        Computes the shortest vector between pos_0 and pos_1 using box 
        dimensions and accounting for periodic boundaries.
        
        Parameters
        ----------
        pos_0 : index or array-like object
            Absolute Cartesian vector position(s) to use as reference point(s).
            If the value can be used as an index, then self.atoms.pos[pos_0]
            is taken.
        pos_1 : index or array-like object
            Absolute Cartesian vector position(s) to find relative to pos_0.
            If the value can be used as an index, then self.atoms.pos[pos_1]
            is taken.

        Returns
        -------
        numpy.ndarray
            The shortest vectors from each pos_0 to pos_1 positions.
        """
        # Test if pos_0 and pos_1 can be used as numpy array indices
        try:
            pos_0 = self.atoms.pos[pos_0]
        except:
            pos_0 = np.asarray(pos_0)
        try:
            pos_1 = self.atoms.pos[pos_1]
        except:
            pos_1 = np.asarray(pos_1)
        
        # Call dvect using self's box and pbc
        vects = dvect(pos_0, pos_1, self.box, self.pbc)
        if len(vects) == 1:
            return vects[0]
        else:
            return vects

    def dmag(self,
             pos_0: Union[int, list, slice, npt.ArrayLike],
             pos_1: Union[int, list, slice, npt.ArrayLike]
             ) -> np.ndarray:
        """
        Computes the shortest distance between pos_0 and pos_1 using box 
        dimensions and accounting for periodic boundaries.
        
        Parameters
        ----------
        pos_0 : index or array-like object
            Absolute Cartesian vector position(s) to use as reference point(s).
            If the value can be used as an index, then self.atoms.pos[pos_0]
            is taken.
        pos_1 : index or array-like object
            Absolute Cartesian vector position(s) to find relative to pos_0.
            If the value can be used as an index, then self.atoms.pos[pos_1]
            is taken.
        
        Returns
        -------
        numpy.ndarray
            The shortest vector magnitude from each pos_0 to pos_1 positions.
        """
        # Test if pos_0 and pos_1 can be used as numpy array indices
        try:
            pos_0 = self.atoms.pos[pos_0]
        except:
            pos_0 = np.asarray(pos_0)
        try:
            pos_1 = self.atoms.pos[pos_1]
        except:
            pos_1 = np.asarray(pos_1)
        
        # Call dvect using self's box and pbc
        vects = dmag(pos_0, pos_1, self.box, self.pbc)
        if len(vects) == 1:
            return vects[0]
        else:
            return vects
    
    def neighborlist(self, **kwargs) -> NeighborList:
        """
        Builds a neighbor list for the system.  The resulting NeighborList
        object is saved to the object as attribute 'neighbors'.
        
        Parameters
        ----------
        cutoff : float, optional
            Radial cutoff distance for identifying neighbors.  Must be given if
            model is not given.
        model : str or file-like object, optional
            Gives the file path or content to load.  If given, no other
            parameters are allowed.
        initialsize : int, optional
            The number of neighbor positions to initially assign to each atom.
            Default value is 20.
        deltasize : int, optional
            Specifies the number of extra neighbor positions to allow each atom
            when the number of neighbors exceeds the underlying array size.
            Default value is 10.
            
        Returns
        -------
        atomman.NeighborList
            The compiled list of neighbors.
        """
        if 'system' in kwargs:
            raise KeyError("Parameter 'system' not allowed")
        else:
            kwargs['system'] = self
        return NeighborList(**kwargs)
    
    def r0(self,
           neighbors: Optional[NeighborList] = None) -> float:
        """
        Identifies the shortest interatomic spacing between atoms in the
        system by comparing the shortest periodic box vector with the
        smallest dmags of all neighbor atoms.
        
        Parameters
        ----------
        neighbors : NeighborList, optional
            A pre-computed NeighborList for the system.  If not given, a new
            NeighborList will be computed using a cutoff distance based on the
            smallest dmag between atom 0 and the rest of the system's atoms.
        
        Returns
        -------
        float
            The shortest interatomic spacing identified.
        """
        
        # Find shortest periodic lattice parameter
        try:
            box_r0 = np.linalg.norm(self.box.vects, axis=1)[self.pbc].min()
        except:
            box_r0 = None
        
        # Find shortest interatomic vector
        if self.natoms > 1:
            
            # Compute neighbors if needed
            if neighbors is None:
                
                # Use r0 for atom 0 as cutoff
                cutoff = self.dmag(0, range(1, self.natoms)).min() * 1.01
                
                # Identify all neighbors
                neighbors = self.neighborlist(cutoff=cutoff)
            
            # Find smallest r0 across all neighbor sets
            atom_r0 = np.inf
            for i in range(self.natoms):
                j = neighbors[i]
                if len(j) > 0:
                    new_r0 = self.dmag(i, j).min()
                    if new_r0 < atom_r0:
                        atom_r0 = new_r0

            if atom_r0 == np.inf:
                atom_r0 = None
        else:
            atom_r0 = None

        if box_r0 is not None and atom_r0 is not None:
            if atom_r0 < box_r0:
                return atom_r0
            else:
                return box_r0
        elif box_r0 is not None:
            return box_r0
        elif atom_r0 is not None:
            return atom_r0
        else:
            raise ValueError('No atoms to compare or periodic boundaries found!')
    
    def supersize(self,
                  a_size: Union[int, Tuple[int, int]],
                  b_size: Union[int, Tuple[int, int]],
                  c_size: Union[int, Tuple[int, int]]) -> System:
        """
        Creates a larger system from a given system by replicating it along the
        system's box vectors.
        
        The multiplier values \\*_size are taken to be integer tuples (m, n) where
        m <= 0 and n >= 0.  The system multiplication works such that if n = -m,
        then the seed system's origin will be at the center of the new system.
        If only one integer is given, then it is assigned to m or n depending on
        its sign, and the other value is taken to be 0.
        
        Parameters
        ----------
        a_size : int or tuple of int
            Single int or two integers specifying replication along the avect
            direction.
        b_size -- int or tuple of int
            Single int or two integers specifying replication along the bvect
            direction.
        c_size -- int or tuple of int
            Single int or two integers specifying replication along the cvect
            direction.
        
        Returns
        -------
        atomman.System
            A new system created by replicating the given seed system according to
            the \\*_size parameters.
          
        """
        # Extract parameters
        sizes = [a_size, b_size, c_size]
        mults = np.array([0, 0, 0], dtype=int)
        vects = self.box.vects
        origin = self.box.origin
        spos = self.atoms_prop('pos', scale=True)
        
        # Check the *_size values
        for i in range(3):
            
            # Change single int to tuple of two int
            if isinstance(sizes[i], (int, np.integer)):
                if sizes[i] > 0:
                    sizes[i] = (0, sizes[i])
                elif sizes[i] < 0:
                    sizes[i] = (sizes[i], 0)
            
            elif isinstance(sizes[i], tuple):
                try:
                    assert len(sizes[i]) == 2, str(len(sizes[i]))
                    assert isinstance(sizes[i][0], (int, np.integer)), str(sizes[i][0])
                    assert sizes[i][0] <= 0, str(sizes[i][0])
                    assert isinstance(sizes[i][1], (int, np.integer)), str(sizes[i][1])
                    assert sizes[i][1] >= 0, str(sizes[i][1])
                except:
                    raise TypeError('Invalid system multipliers')
            else:
                raise TypeError('Invalid system multipliers')
            
            # Calculate full multipliers 
            mults[i] = sizes[i][1] - sizes[i][0]
            if mults[i] == 0:
                raise ValueError('Cannot multiply system dimension by zero')
                
            # Scale box and first set of positions accordingly
            spos[:,i] /= mults[i]
            origin += vects[i] * sizes[i][0]
            vects[i] *= mults[i]
            
        # Initialize new Box and Atoms
        box = Box(vects=vects, origin=origin)
        natoms = self.natoms * mults[0] * mults[1] * mults[2]
        atoms = Atoms(natoms=natoms)
        
        # Copy over all property values (except pos) using numpy broadcasting
        for key in self.atoms_prop():
            if key == 'pos':
                continue
        
            # Get old array
            old = self.atoms.view[key]
            
            # Create new array and broadcast old to it
            new = np.empty((mults[0] * mults[1] * mults[2],) + old.shape, dtype = old.dtype)
            new[:] = old
            
            # Reshape new and save to atoms
            new_shape = new.shape
            new_shape = (new_shape[0] * new_shape[1],) + new_shape[2:] # pylint: disable=unsubscriptable-object
            atoms.view[key] = np.array(new.reshape(new_shape))
        
        # Expand spos using broadcasting
        new_spos = np.empty((mults[0] * mults[1] * mults[2],) + spos.shape)
        new_spos[:] = spos
        
        # Reshape spos
        new_shape = new_spos.shape
        new_shape = (new_shape[0]*new_shape[1],) + new_shape[2:] # pylint: disable=unsubscriptable-object
        new_spos = new_spos.reshape(new_shape)
        
        # Use broadcasting to create arrays to add to spos
        test = np.empty(mults[0] * self.natoms)
        test.shape = (self.natoms, mults[0])
        test[:] = np.arange(mults[0])
        x = test.T.flatten()

        test = np.empty(mults[1] * len(x))
        test.shape = (len(x), mults[1])
        test[:] = np.arange(mults[1])
        y = test.T.flatten()
        test.shape = (mults[1], len(x))
        test[:] = x
        x = test.flatten()

        test = np.empty(mults[2] * len(x))
        test.shape = (len(x), mults[2])
        test[:] = np.arange(mults[2])
        z = test.T.flatten()
        test.shape = (mults[2], len(x))
        test[:] = x
        x = test.flatten()
        test[:] = y
        y = test.flatten()
        
        # xyz is displacement values to add to spos
        xyz = (np.hstack((x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis])) 
             * np.array([1 / mults[0], 1 / mults[1], 1 / mults[2]]))
        
        # Save pos values, return new System
        atoms.view['pos'] = new_spos + xyz
        
        return System(box=box, atoms=atoms, scale=True, symbols=self.symbols)
    
    def rotate(self,
               uvws: npt.ArrayLike,
               tol: Union[float, list, None] = None,
               return_transform: bool = False) -> System:
        """
        Transforms a System representing a periodic crystal cell from a standard
        orientation to a specified orientation. Note: if hexagonal indices are
        given, the vectors will be reduced to the smallest uvw integer
        representation.
        
        Parameters
        ----------
        uvws : array-like object
            A (3, 3) array of the Miller crystal vectors or a (3, 4) array of
            Miller-Bravais hexagonal crystal vectors to use in transforming the
            system.  Values must be integers.
        tol : list or float, optional
            Tolerance parameter used in determining which atoms are inside the
            box.  Multiple values can be given as the identification may
            occasionally fail for a given ucell and tol.  Default behavior will
            try tol values ranging from 1e-4 to 1e-8.
        return_transform : bool, optional
            Indicates if the transformation matrix associated with the
            rotation is returned.  Default value is False.
        
        Returns
        -------
        atomman.System
            A new fully periodic system rotated and transformed according to the
            uvws crystal vectors.
        transform : np.ndarray
            The transformation matrix associated with the rotation.
            Returned if return_transform is True.
        """
        
        if tol is None:
            tol = [1e-4, 1e-5, 1e-6, 1e-7]
        else:
            tol = aslist(tol)
        
        uvws = np.asarray(uvws)
        
        # Convert uvws from Miller-Bravais to Miller indices if needed
        if uvws.shape == (3, 4):
            if self.box.ishexagonal():
                uvws = miller.vector4to3(uvws)
            else:
                raise ValueError('hexagonal indices only work on hexagonal systems')
        
        # Check uvws shape and values
        if uvws.shape != (3, 3):
            raise ValueError('Invalid uvws crystal indices shape')

        int_uvws = np.asarray(np.rint(uvws), dtype='int64')
        if np.allclose(uvws, int_uvws):
            uvws = int_uvws
        else:
            raise ValueError('Rotation uvws must be integer values')
        
        # No rotation shortcut
        if np.all(uvws == np.eye(3, dtype='int64')):
            newsystem = deepcopy(self)
        
        else:
            # Get natoms and volume of system
            natoms = self.natoms
            volume = self.box.volume
            
            # Convert uvws to Cartesian units and compute new volume and natoms
            newvects = miller.vector_crystal_to_cartesian(uvws, box=self.box)
            newvolume = np.abs(newvects[0].dot(np.cross(newvects[1], newvects[2])))
            newnatoms = int(round(newvolume / volume) * natoms)
            
            # Check new values
            if newnatoms == 0:
                raise ValueError('New box has no atoms/volume: vectors are parallel or planar')
            
            # Identify box corners of new system wrt uvws
            corners = np.empty((8,3), dtype='int64')
            corners[0] = np.zeros(3)
            corners[1] = uvws[0]
            corners[2] = uvws[1]
            corners[3] = uvws[2]
            corners[4] = uvws[0] + uvws[1]
            corners[5] = uvws[0] + uvws[2]
            corners[6] = uvws[1] + uvws[2]
            corners[7] = uvws[0] + uvws[1] + uvws[2]
            
            # Create a supercell of system that contains all box corners
            a_mults = (corners[:,0].min()-1, corners[:,0].max()+1)
            b_mults = (corners[:,1].min()-1, corners[:,1].max()+1)
            c_mults = (corners[:,2].min()-1, corners[:,2].max()+1)
            system2 = self.supersize(a_mults, b_mults, c_mults)
            
            # Change system.box.vects to newvects
            system2.box_set(vects=newvects, scale=False)
            
            search_success = False
            for atol in tol:
                
                spos = system2.atoms_prop('pos', scale=True)
                
                # Round atom positions near box boundaries to the boundaries
                spos[np.isclose(spos, 0.0, atol=atol)] = 0.0
                spos[np.isclose(spos, 1.0, atol=atol)] = 1.0
            
                # Identify all atoms whose positions are 0 <= x < 1
                aindex = np.where(((spos[:, 0] >= 0.0) & (spos[:, 0] < 1.0)
                                 & (spos[:, 1] >= 0.0) & (spos[:, 1] < 1.0)
                                 & (spos[:, 2] >= 0.0) & (spos[:, 2] < 1.0)))
                
                # Check if number of atoms identified matches the expected number
                if len(aindex[0]) == newnatoms:
                    search_success = True
                    break
            
            if not search_success:
                raise ValueError(f'Filtering failed: {newnatoms} atoms expected, {len(aindex[0])} found')
            
            # Make newsystem by cutting out all atoms in system2 outside boundaries
            newsystem = System(atoms=system2.atoms[aindex], box=system2.box, symbols=self.symbols)
        
        # Return normalized system
        return newsystem.normalize(return_transform=return_transform)
    
    def normalize(self,
                  style: str = 'lammps',
                  return_transform: bool = False
                  ) -> Union[System, Tuple[System, np.ndarray]]:
        """
        Normalizes a system's box vectors and atom positions to be compatible
        with simulation codes.
        
        Parameters
        ----------
        style : str, optional
            Indicates the normalization style to use.  Default (and only
            current option) is 'lammps'.
        return_transform : bool, optional
            Indicates if the transformation matrix associated with the
            normalization is returned.  Default value is False.
        Returns
        -------
        newsystem : atomman.System
            A new system that has been normalized.
        transform : np.ndarray
            The transformation matrix associated with the normalization.
            Returned if return_transform is True.
        """
        if style == 'lammps':
            return lmp_normalize(self, return_transform=return_transform)
        else:
            raise ValueError("Unknown style (only 'lammps' is currently supported)")
    
    def dump(self, style: str, **kwargs) -> Optional[Any]:
        """
        Convert a System to another format.
    
        Parameters
        ----------
        style : str
            Indicates the format of the content to dump the atomman.System as.
        kwargs
            Any extra keyword arguments to pass to the underlying dump methods.
            
        Returns
        -------
        str, object or tuple
            Any content returned by the underlying dump methods.
        """
        return dump(style, self, **kwargs)

    def model(self,
              box_unit: Optional[str] = None,
              prop_name: Optional[list] = None,
              unit: Optional[list] = None,
              prop_unit: Optional[dict] = None) -> DM:
        """
        Generates a data model for the System object.

        Parameters
        ----------
        box_unit : str, optional
            Length unit to use for the box. Default value is 'angstrom'.
        prop_name : list, optional
            The Atoms properties to include.  If neither prop_name nor prop_unit
            are given, all system properties will be included.
        unit : list, optional
            Lists the units for each prop_name as stored in the table.  For a
            value of None, no conversion will be performed for that property.  For
            a value of 'scaled', the corresponding table values will be taken in
            box-scaled units.  If neither unit nor prop_units given, pos will be
            given in Angstroms and all other values will not be converted.
        prop_unit : dict, optional
            dictionary where the keys are the property keys to include, and
            the values are units to use. If neither unit nor prop_units given, 
            pos will be given in Angstroms and all other values will not be
            converted.

        Returns
        -------
        DataModelDict.DataModelDict
            A JSON/XML data model for the current System object. 
        """

        # Initialize DataModelDict
        model = DM()
        model['atomic-system'] = DM()

        # Add box
        model['atomic-system']['box'] = self.box.model(length_unit=box_unit)['box']
        
        # Add pbc
        model['atomic-system']['periodic-boundary-condition'] = self.pbc.tolist()
        
        # Add symbols
        for symbol in self.symbols:
            model['atomic-system'].append('atom-type-symbol', symbol)

        # Add masses
        addmasses = False
        for mass in self.masses:
            if mass is not None:
                addmasses = True
                break
        if addmasses:
            for mass in self.masses:
                model['atomic-system'].append('atom-type-mass', mass)
        
        # Add atoms
        model['atomic-system']['atoms'] = amodel = self.atoms.model(prop_name=prop_name,
                                                                    unit=unit, 
                                                                    prop_unit=prop_unit)['atoms']
        
        # Scale properties if needed
        for prop in amodel.aslist('property'):
            if prop['data'].get('unit', None) == 'scaled':
                prop['data'] = uc.model(self.box.position_cartesian_to_relative(uc.value_unit(prop['data'])), units='scaled') 
        
        return model
