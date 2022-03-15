# coding: utf-8
# Standard Python libraries
from __future__ import annotations
import io
from copy import deepcopy
from collections import OrderedDict
from typing import Any, Optional, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# https://pandas.pydata.org/
import pandas as pd

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
import atomman.unitconvert as uc
from ..tools import indexstr

class Atoms(object):
    """
    Class for representing a collection of atoms.
    """
    
    class PropertyDict(OrderedDict):
        """Extends OrderedDict to work with Atoms"""
        
        def __init__(self, host: Atoms):
            """PropertyDict needs to know its host Atoms object."""
            self.__host = host
            super(Atoms.PropertyDict, self).__init__()
        
        def __setitem__(self, key: str, value: np.ndarray):
            """
            Modifies OrderedDict__setitem__() such that
            1. All values are converted to numpy.ndarrays.
            2. For new keys, checks are done that data is compatible with host.natoms.
            3. New keys are also added as attributes to host Atoms object, if allowed.
            4. For existing keys, new values are saved over the old ones as opposed to
            only changing the name assignment.

            Parameters
            ----------
            key : str
                The property key name to assign values to.
            value : numpy.ndarray
                The values to assign.
            """
            # Shortcut to host
            host = self.__host
            
            # Python 2: change to unicode if needed.
            try:
                key = key.decode('UTF-8')
            except:
                pass
            
            # Convert to numpy.ndarray if needed
            value = np.asarray(value)
            
            # Broadcast if needed and allowed
            if value.shape == ():
                value = np.array(np.broadcast_to(value, (host.natoms,) + value.shape))
            elif value.shape[0] == 1:
                value = np.array(np.broadcast_to(value, (host.natoms,) + value.shape[1:]))
            elif value.shape[0] != host.natoms:
                raise ValueError('First dimension of value must be 1 or natoms')
            
            # Check that atype values are 1 or greater
            if key == 'atype' and len(value) > 0 and np.min(value) < 1:
                raise ValueError('atype values must be >= 1')

            # If key is already assigned, save value over existing values
            if key in self.keys():
                self[key][:] = value
            
            # Otherwise, set new item and try assigning attribute to host
            else:
                super(Atoms.PropertyDict, self).__setitem__(key, value)
                try:
                    assert key not in dir(host)
                    super(Atoms, host).__setattr__(key, value) # pylint: disable=bad-super-call
                except:
                    pass
    
    def __init__(self,
                 natoms: Optional[int] = None,
                 atype: Union[int, npt.ArrayLike, None] = None,
                 pos: Optional[npt.ArrayLike] = None,
                 prop: Optional[dict] = None,
                 model: Union[str, io.IOBase, DM, None] = None, 
                 safecopy: bool = False,
                 **kwargs):
        """
        Class initializer.
        
        Parameters
        ==========
        natoms : int, optional
            The number of atoms.  If not given, will be inferred from other
            parameters if possible, or set to 1 if not possible.
        atype : int or list/ndarray of int, optional
            The integer atomic types to assign to all atoms.  Default is to
            set all atypes to 1.
        pos : list/ndarray of float, optional
            The atomic positions to assign to all atoms.  Default is to set
            each atom's position to [0,0,0].
        model : str or DataModelDict, optional
            File path or content of a JSON/XML data model containing all
            atom information.  Cannot be given with any other parameters.
        prop : dict, optional
            Dictionary containing all per-atom properties.  Can be used
            instead of atype, pos, and kwargs.  This is for backwards
            compatibility support with atomman version 1.
        safecopy : bool, optional
            Flag indicating if values are to be copied before setting.  For
            property values given as numpy arrays, direct setting (False,
            default) may result in the Atoms' property pointing to the original
            numpy array.  Using safecopy=True deep copies the property values
            before setting to avoid this.  Note that safecopy=True may be
            considerably slower for large numbers of atoms and/or properties.
        kwargs : any
            All additional key/value pairs are assigned as properties.
        """
        
        # Check for model
        if model is not None:
            try:
                assert natoms is None
                assert atype is None
                assert pos is None
                assert prop is None
                assert len(kwargs) == 0
            except:
                raise ValueError('model cannot be given with any other parameters')
            
            # Extract natoms and properties from data model
            model = DM(model).find('atoms')
            natoms = model['natoms']
            prop = OrderedDict()
            for propmodel in model.aslist('property'):
                prop[propmodel['name']] = uc.value_unit(propmodel['data'])
        
        # Check for prop dictionary
        if prop is not None:
            if atype is not None and pos is not None and len(kwargs) > 0:
                raise ValueError('prop dict cannot be given with keyword properties')
            
            # Divide prop into keyword arguments
            atype = prop.pop('atype', None)
            pos = prop.pop('pos', None)
            kwargs = prop
        
        # Check atype parameter values
        if atype is not None:
            atype = np.asarray(atype)
            
            # Handle single atype
            if atype.ndim == 0:
                natoms_atype = 1
            
            # Handle list of atype
            elif atype.ndim == 1:
                natoms_atype = atype.shape[0]
            else:
                raise ValueError('Only one atype per atom allowed')
        else:
            atype = np.array([1], dtype='uint64')
            natoms_atype = 1
        
        # Check pos parameter values
        if pos is not None:
            pos = np.asarray(pos)
            
            # Handle single pos
            if pos.ndim == 1:
                if pos.shape[0] == 3:
                    natoms_pos = 1
                else:
                    raise ValueError('pos per atom must be 3-dimensional')
            
            # Handle list of pos
            elif pos.ndim == 2:
                if pos.shape[1] == 3:
                    natoms_pos = pos.shape[0]
                else:
                    raise ValueError('pos per atom must be 3-dimensional')
            else:
                raise ValueError('too many dimensions for pos')
        else:
            pos = np.zeros((1, 3), dtype='float64')
            natoms_pos = 1
        
        # Check natoms parameter values
        if natoms is not None:
            natoms = int(natoms)
            if ((natoms_atype == 1 or natoms_atype == natoms)
                and (natoms_pos == 1 or natoms_pos == natoms)):
                pass
            else:
                raise ValueError('natoms and length of atype/pos not compatible')
        else:
            if natoms_atype == natoms_pos:
                natoms = natoms_atype
            elif natoms_atype == 1:
                natoms = natoms_pos
            elif natoms_pos == 1:
                natoms = natoms_atype
            else:
                raise ValueError('lengths of atype and pos not compatible')
            
        # Initialize underlying private class attributes
        super(Atoms, self).__setattr__('_Atoms__natoms', natoms)
        super(Atoms, self).__setattr__('_Atoms__view', Atoms.PropertyDict(self))
        super(Atoms, self).__setattr__('_Atoms__dir', deepcopy(dir(self)))
        
        # Set properties
        if safecopy:
            self.view['atype'] = deepcopy(atype)
            self.view['pos'] = deepcopy(pos)
            for key, value in kwargs.items():
                self.view[key] = deepcopy(value)
        else:
            self.view['atype'] = atype
            self.view['pos'] = pos
            for key, value in kwargs.items():
                self.view[key] = value
    
    def model(self,
              prop_name: Optional[list] = None,
              unit: Optional[list] = None,
              prop_unit: Optional[dict] = None) -> DM:
        """
        Generates a data model for the Atoms object.
        
        Parameters
        ----------
        prop_name : list, optional
            The Atoms properties to include.  If neither prop_name nor prop_unit
            are given, all system properties will be included.
        unit : list, optional
            Lists the units for each prop_name as stored in the table.  For a
            value of None, no conversion will be performed for that property.  
            If neither unit nor prop_units given, pos will be
            given in Angstroms and all other values will not be converted.
        prop_unit : dict, optional
            dictionary where the keys are the property keys to include, and
            the values are units to use. If neither unit nor prop_units given, 
            pos will be given in Angstroms and all other values will not be
            converted.

        Returns
        -------
        DataModelDict.DataModelDict
            A JSON/XML data model for the current Atoms object.  
        """
        
        # Set prop_unit if needed
        if prop_unit is None:
            if prop_name is None:
                prop_name = self.prop()
            if unit is None:
                unit = [None for i in range(len(prop_name))]
            
            if len(unit) != len(prop_name):
                raise ValueError('')
            prop_unit = {}
            for p, u  in zip(prop_name, unit):
                prop_unit[p] = u
        
        elif prop_name is not None or unit is not None:
            raise ValueError('prop_unit cannot be given with prop_name or unit')
        
        # Set default pos unit
        if 'pos' in prop_unit and prop_unit['pos'] is None:
            prop_unit['pos'] = 'angstrom'
        
        # Generate DataModelDict
        model = DM()
        model['atoms'] = DM()
        model['atoms']['natoms'] = self.natoms
        for prop in prop_unit:
            unit = prop_unit.get(prop, None)
            propmodel = DM()
            propmodel['name'] = prop
            propmodel['data'] = uc.model(self.prop(prop), unit)
            model['atoms'].append('property', propmodel)
        
        return model

    def __str__(self) -> str:
        """string output of Atoms data"""
        atype = self.atype # pylint: disable=no-member
        pos = self.pos # pylint: disable=no-member
        lines = ['per-atom properties = ' + str(self.prop())]
        lines.append('     id |   atype |  pos[0] |  pos[1] |  pos[2]')
        for i in range(self.natoms):
            lines.append('%7i | %7i | %7.3f | %7.3f | %7.3f' % (i, atype[i], pos[i,0], pos[i,1], pos[i,2]))
        
        return '\n'.join(lines)
    
    def __setattr__(self, name: str, value: Any):
        """Control setting of user-defined attributes"""
        if not hasattr(self, name) or name in self.view:
            self.view[name] = value
        else:
            super(Atoms, self).__setattr__(name, value)
    
    def __intslice(self, intnum):
        if intnum == -1:
            return slice(intnum, None)
        else:
            return slice(intnum, intnum+1)
    
    def __deepcopy__(self, memo) -> Atoms:
        """Properly handle deepcopy"""
        d = OrderedDict()
        atype = deepcopy(self.view['atype'])
        pos = deepcopy(self.view['pos'])
        for key in self.view:
            if key not in ['atype', 'pos']:
                d[key] = deepcopy(self.view[key])
        return Atoms(atype=atype, pos=pos, **d)
    
    def __getitem__(self, index) -> Atoms:
        """Index getting of Atoms."""
        view = OrderedDict()
        if isinstance(index, (int, np.integer)):
            index = self.__intslice(index)
        for key in self.view.keys():
            view[key] = self.view[key][index]
        return Atoms(**view)
        
    def __setitem__(self, index, value):
        """Index setting of Atoms."""
        try:
            assert isinstance(value, Atoms)
            assert sorted(value.view.keys()) == sorted(self.view.keys())
        except:
            raise ValueError('Can only set Atoms with matching properties')
        
        if isinstance(index, (int, np.integer)):
            index = self.__intslice(index)
            
        for key in self.view.keys():
            self.view[key][index] = value.view[key]
    
    def __len__(self) -> int:
        """len of atoms = natoms."""
        return self.natoms
    
    @property
    def natoms(self) -> int:
        """int : The number of atoms in the Atoms class."""
        return self.__natoms 
    
    @property
    def atypes(self) -> tuple:
        """tuple : List of int atom types."""
        return tuple(range(1, self.natypes+1))
    
    @property
    def natypes(self) -> int:
        """int : The number of atom types in the Atoms class."""
        if np.min(self.atype) < 1: 
            raise ValueError('atype values < 1 not allowed')

        return int(np.max(self.atype)) 
        
    @property
    def view(self) -> PropertyDict:
        """PropertyDict : All assigned per-atom properties."""
        return self.__view 
    
    def prop(self,
             key: Optional[str] = None,
             index: Union[int, list, slice, None] = None,
             value: Optional[Any] = None,
             a_id: Optional[int] = None
             ) -> Union[list, Atoms, np.ndarray, None]:
        """
        Accesses the per-atom properties for controlled getting and setting.
        For getting values, prop() always returns a copy of the underlying
        data as opposed to the data itself.
        
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
        
        raises
        ------
        ValueError
            If a_id and index are both given, or only value is given.
        """
        
        # Handle a_id
        if a_id is not None:
            if index is not None:
                raise ValueError('a_id and index cannot both be given')
            index = a_id
        
        # Get values if value is None
        if value is None:
            if key is None:
                # Return list of property names when no parameters are given
                if index is None:
                    return list(self.view.keys())
                # Return copy of slice if only index is given
                else:
                    return deepcopy(self[index])
            # If key is given, return copy of property values
            else:
                if index is None:
                    return deepcopy(self.view[key])
                else:
                    return deepcopy(self.view[key][index])
        
        # Set values if value is given
        else:
            # If no key, set value to atoms (or a slice)
            if key is None:
                if not isinstance(value, Atoms):
                    raise TypeError('If key is None, value must be instance of atomman.Atoms')
                if index is None:
                    self[:] = value
                else:
                    self[index] = value
            
            # If key is given, set copy of value to property
            else:
                if index is None:
                    self.view[key] = deepcopy(value)
                else:
                    self.view[key][index] = value
    
    def prop_atype(self,
                   key: str,
                   value: Any,
                   atype: Optional[int] = None):
        """
        Allows for per-atom properties to be assigned according to
        Atoms.atypes.

        Parameters
        ----------
        key : str
            Per-atom property name.
        value : list, any
            Property value(s) to assign.  If atype is not given, this should be
            an object of length Atoms.natypes. Otherwise, should be a single per-atom
            value.
        atype : int, optional
            A specific atype to assign value to.

        Raises
        ------
        ValueError
            If length of value does not match Atoms.natypes or atype is not in
            Atoms.atypes.
        """

        # Set values across all atype values
        if atype is None:
            value = np.asarray(value)
            if len(value) >= self.natypes:
                self.view[key] = value[self.atype - 1] # pylint: disable=no-member
            else:
                raise ValueError('length of value less than natypes')
        
        # Set values for only one atype
        else:
            if atype in self.atypes:
                if key not in self.prop():
                    self.view[key] = np.zeros_like(value)
                self.view[key][self.atype==atype] = value # pylint: disable=no-member
            else:
                raise ValueError('atype not found')
    
    def df(self) -> pd.DataFrame:
        """
        Returns a pandas.DataFrame of all atomic properties.  Multi-dimensional
        per-atom data will be converted into multiple table columns.
        """
        
        # Initialize new dictionary of values
        values = OrderedDict()
        for key in self.view.keys():
            
            # Flatten multidimensional arrays
            for index, istr in indexstr(self.view[key].shape[1:]):
                newkey = key + istr
                
                # Copy values over
                if index == ():
                    values[newkey] = self.view[key]
                else:
                    values[newkey] = self.view[key][(Ellipsis, ) + index]
        
        # Return DataFrame
        return pd.DataFrame(values)
    
    def extend(self, value: Union[Atoms, int]) -> Atoms:
        """
        Allows additional atoms to be added to the end of the atoms list.
        
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
        
        Returns
        -------
        atomman.Atoms
            A new Atoms object containing all atoms and properties of the
            current object plus the additional atoms.
        """
        
        # Handle different value types
        if isinstance(value, (int, np.integer)):
            natoms = value
            atoms = Atoms(natoms=natoms)
        elif isinstance(value, Atoms):
            natoms = value.natoms
            atoms = value        
        else:
            raise TypeError('can only add Atoms or an int # of atoms')
        
        # Generate newatoms by copying values of self + natoms extras
        index = list(range(self.natoms)) + [0 for i in range(natoms)]
        newatoms = self[index]
        
        # Create empty values for atoms.props not in newatoms (self)
        for prop in atoms.prop():
            if prop not in newatoms.prop():
                newatoms.view[prop] = np.zeros((newatoms.natoms, ) + atoms.view[prop][0].shape, dtype=atoms.view[prop][0].dtype)
        
        # Copy values to the extra atoms in newatoms
        for prop in newatoms.prop():
            if prop in atoms.prop():
                newatoms.view[prop][self.natoms:] = atoms.view[prop]
            else:
                newatoms.view[prop][self.natoms:] = np.zeros((natoms, ) + self.view[prop][0].shape, dtype=self.view[prop][0].dtype)
        
        return newatoms