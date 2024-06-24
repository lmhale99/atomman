# coding: utf-8

# Standard Python imports
import io
from copy import deepcopy
from typing import Optional, Union, Tuple

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_query, load_value

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

class Dislocation(Record):
    """
    Class for representing dislocation records, which collect the parameters
    necessary for atomman to generate a particular dislocation type.
    """
    def __init__(self,
                 model: Union[str, io.IOBase, DM, None] = None,
                 name: Optional[str] = None,
                 database = None):
        """
        Initializes a Record object for a given style.
        
        Parameters
        ----------
        model : str, file-like object, DataModelDict
            The contents of the record.
        name : str, optional
            The unique name to assign to the record.  If model is a file
            path, then the default record name is the file name without
            extension.
        database : yabadaba.Database, optional
            Allows for a default database to be associated with the record.
        """
        super().__init__(model=model, name=name, database=database)

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'dislocation'

    @property
    def xsd_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsd schema"""
        return ('atomman.library.xsd', f'{self.style}.xsd')

    @property
    def xsl_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsl transformer"""
        return ('atomman.library.xsl', f'{self.style}.xsl')

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'dislocation'

    ####################### Define Values and attributes #######################

    def _init_value_objects(self) -> list:
        """
        Method that defines the value objects for the Record.  This should
        1. Call the method's super() to get default Value objects.
        2. Use yabadaba.load_value() to build Value objects that are set to
           private attributes of self.
        3. Append the list returned by the super() with the new Value objects.

        Returns
        -------
        value_objects: A list of all value objects.
        """
        value_objects = super()._init_value_objects()
        
        self.__key = load_value('str', 'key', self, valuerequired=True)
        self.__id = load_value('str', 'id', self, valuerequired=True)
        self.__url = load_value('str', 'url', self, modelpath='URL')
        self.__family = load_value('str', 'family', self, valuerequired=True,
                                   modelpath='system-family')
        self.__family_url = load_value('str', 'family_url', self,
                                   modelpath='system-family-URL')
        self.__character = load_value('str', 'character', self, valuerequired=True,
                                      allowedvalues=['screw', 'edge', 'mixed'])
        self.__slip_hkl = load_value('miller', 'slip_hkl', self, valuerequired=True,
                                     modelpath='calculation-parameter.slip_hkl',
                                     bracket='()')
        self.__line_uvw = load_value('miller', 'line_uvw', self, valuerequired=True,
                                     modelpath='calculation-parameter.line_uvw',
                                     bracket='[]')
        self.__burgers_uvw = load_value('miller', 'burgers_uvw', self, valuerequired=True,
                                        modelpath='calculation-parameter.burgers_uvw',
                                        bracket='[]')
        self.__m = load_value('unitvector', 'm', self,
                                  modelpath='calculation-parameter.m')
        self.__n = load_value('unitvector', 'n', self,
                                  modelpath='calculation-parameter.n')
        self.__shift = load_value('vector', 'shift', self,
                                  modelpath='calculation-parameter.shift')
        self.__shiftscale = load_value('bool', 'shiftscale', self,
                                       modelpath='calculation-parameter.shiftscale')
        self.__shiftindex = load_value('int', 'shiftindex', self,
                                       modelpath='calculation-parameter.shiftindex')

        value_objects.extend([
            self.__key, self.__id, self.__url, self.__family, self.__family_url, self.__character,
            self.__slip_hkl, self.__line_uvw, self.__burgers_uvw, self.__m, self.__n, 
            self.__shift, self.__shiftscale, self.__shiftindex])

        return value_objects

    @property
    def key(self) -> str:
        """str : A UUID4 key assigned to the record"""
        return self.__key.value

    @key.setter
    def key(self, val: str):
        self.__key.value = val

    @property
    def id(self) -> str:
        """str : A unique id assigned to the record"""
        return self.__id.value

    @id.setter
    def id(self, val: str):
        self.__id.value = val

    @property
    def url(self) -> Optional[str]:
        """str : A URL where a copy of the record can be found"""
        return self.__url.value

    @url.setter
    def url(self, val: Optional[str]):
        self.__url.value = val

    @property
    def family(self) -> str:
        """str : The prototype/reference id the defect is defined for"""
        return self.__family.value

    @family.setter
    def family(self, val: str):
        self.__family.value = val

    @property
    def family_url(self) -> Optional[str]:
        """str : A URL where a copy of the family record can be found"""
        return self.__family_url.value

    @family_url.setter
    def family_url(self, val: Optional[str]):
        self.__family_url.value = val

    @property
    def character(self) -> str:
        """str : The dislocation's character"""
        return self.__character.value

    @character.setter
    def character(self, val: str):
        self.__character.value = val

    @property
    def slip_hkl(self) -> np.ndarray:
        """numpy.ndarray : The dislocation's slip plane"""
        return self.__slip_hkl.value

    @slip_hkl.setter
    def slip_hkl(self, val: Union[str, npt.ArrayLike]):
        self.__slip_hkl.value = val
    
    @property
    def slip_hkl_str(self) -> str:
        """str: The string representation of slip_hkl"""
        return self.__slip_hkl.str
    
    @property
    def line_uvw(self) -> np.ndarray:
        """numpy.ndarray : The dislocation's slip plane"""
        return self.__line_uvw.value

    @line_uvw.setter
    def line_uvw(self, val: Union[str, npt.ArrayLike]):
        self.__line_uvw.value = val

    @property
    def line_uvw_str(self) -> str:
        """str: The string representation of line_uvw"""
        return self.__line_uvw.str

    @property
    def burgers_uvw(self) -> np.ndarray:
        """numpy.ndarray : The dislocation's Burgers vector"""
        return self.__burgers_uvw.value

    @burgers_uvw.setter
    def burgers_uvw(self, val: Union[str, npt.ArrayLike]):
        self.__burgers_uvw.value = val

    @property
    def burgers_uvw_str(self) -> str:
        """str: The string representation of burgers_uvw"""
        return self.__burgers_uvw.str
    
    @property
    def m(self) -> np.ndarray:
        """np.ndarray: The Cartesian unit vector to align with the dislocation solution's m-axis"""
        return self.__m.value

    @m.setter
    def m(self, val: Union[str, npt.ArrayLike]):
        self.__m.value = val

    @property
    def n(self) -> np.ndarray:
        """np.ndarray: The Cartesian unit vector to align with the dislocation solution's n-axis"""
        return self.__n.value

    @n.setter
    def n(self, val: Union[str, npt.ArrayLike]):
        self.__n.value = val

    @property
    def shift(self) -> Optional[np.ndarray]:
        """np.ndarray or None: A rigid body shift to apply to the atoms in the system before inserting the dislocation"""
        return self.__shift.value
    
    @shift.setter
    def shift(self, val: Union[str, npt.ArrayLike, None]):
        self.__shift.value = val

    @property
    def shiftscale(self) -> Optional[bool]:
        """bool or None: Indicates if shift is absolute Cartesian (False) or to be scaled relative to the rcell (True)"""
        return self.__shiftscale.value
    
    @shiftscale.setter
    def shiftscale(self, val: Optional[bool]):
        self.__shiftscale.value = val

    @property
    def shiftindex(self) -> Optional[int]:
        """int or None: The shift index to use for positioning the slip plane between planes of atoms"""
        return self.__shiftindex.value

    @shiftindex.setter
    def shiftindex(self, val: Optional[int]):
        self.__shiftindex.value = val

    @property
    def parameters(self) -> dict:
        """dict : Defect parameters for atomman structure generator"""
        p = {}
        p['slip_hkl'] = self.slip_hkl_str
        p['line_uvw'] = self.line_uvw_str
        p['burgers_uvw'] = self.burgers_uvw_str
        p['m'] = f"{self.m[0]} {self.m[1]} {self.m[2]}"
        p['n'] = f"{self.n[0]} {self.n[1]} {self.n[2]}"
        if self.shift is not None:
            p['shift'] = f"{self.shift[0]} {self.shift[1]} {self.shift[2]}"
        if self.shiftscale is not None:
            p['shiftscale'] = str(self.shiftscale)
        if self.shiftindex is not None:
            p['shiftindex'] = str(self.shiftindex)

        return p
