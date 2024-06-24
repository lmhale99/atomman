# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Union, Tuple

import numpy as np
import numpy.typing as npt

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_value

class FreeSurface(Record):
    """
    Class for representing free_surface records, which collect the parameters
    necessary for atomman to generate a particular free surface.
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
        return 'free_surface'

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
        return 'free-surface'

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
        self.__hkl = load_value('miller', 'hkl', self, valuerequired=True,
                                modelpath='calculation-parameter.hkl',
                                bracket='()')
        self.__shiftindex = load_value('int', 'shiftindex', self,
                                       modelpath='calculation-parameter.shiftindex')
        self.__cutboxvector = load_value('str', 'cutboxvector', self, valuerequired=True,
                                         modelpath='calculation-parameter.cutboxvector',
                                         defaultvalue='c', allowedvalues=['a', 'b', 'c'])

        value_objects.extend([
            self.__key, self.__id, self.__url, self.__family, self.__family_url,
            self.__hkl, self.__shiftindex, self.__cutboxvector])

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
    def family_url(self) -> Optional[str]:
        """str : A URL where a copy of the family record can be found"""
        return self.__family_url.value

    @family_url.setter
    def family_url(self, val: Optional[str]):
        self.__family_url.value = val

    @property
    def family(self) -> str:
        """str : The prototype/reference id the defect is defined for"""
        return self.__family.value

    @family.setter
    def family(self, val: str):
        self.__family.value = val

    @property
    def hkl(self) -> np.ndarray:
        """numpy.ndarray : The Miller-(Bravais) surface plane"""
        return self.__hkl.value
    
    @hkl.setter
    def hkl(self, val: Union[str, npt.ArrayLike]):
        self.__hkl.value = val

    @property
    def hkl_str(self) -> str:
        """str: The string representation of hkl"""
        return self.__hkl.str

    @property
    def shiftindex(self) -> int:
        """int: The shift index to use for positioning the surface plane between planes of atoms"""
        return self.__shiftindex.value
    
    @shiftindex.setter
    def shiftindex(self, val: Optional[int]):
        self.__shiftindex.value = val

    @property
    def cutboxvector(self) -> str:
        """str: The box vector (a, b, or c) that the plane is cutting along"""
        return self.__cutboxvector.value

    @cutboxvector.setter
    def cutboxvector(self, val: str):
        self.__cutboxvector.value = val

    @property
    def parameters(self) -> dict:
        """dict : Defect parameters for atomman structure generator"""
        p = {}
        p['hkl'] = self.hkl_str
        if self.shiftindex is not None:
            p['shiftindex'] = str(self.shiftindex)
        p['cutboxvector'] = self.cutboxvector

        return p