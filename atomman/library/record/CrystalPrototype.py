# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Union, Tuple

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_query, load_value

# atomman imports
from ... import System

class CrystalPrototype(Record):
    """
    Class for representing crystal_prototype records that describe common
    crystal prototypes.
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
        # Init properties
        self.__ucell = None

        # Call super init
        super().__init__(model=model, name=name, database=database)
    
    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'crystal_prototype'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'crystal-prototype'

    @property
    def xsd_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsd schema"""
        return ('atomman.library.xsd', f'{self.style}.xsd')

    @property
    def xsl_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsl transformer"""
        return ('atomman.library.xsl', f'{self.style}.xsl')

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
        self.__commonname = load_value('str', 'commonname', self, valuerequired=True)
        self.__prototype = load_value('str', 'prototype', self, valuerequired=True)
        self.__pearson = load_value('str', 'pearson', self, valuerequired=True,
                                    modelpath='Pearson-symbol')
        self.__strukturbericht = load_value('str', 'strukturbericht', self, valuerequired=True,
                                            modelpath='Strukturbericht')
        self.__sg_number = load_value('int', 'sg_number', self, valuerequired=True,
                                      modelpath='space-group.number')
        self.__sg_hm = load_value('str', 'sg_hm', self, valuerequired=True,
                                  modelpath='space-group.Hermann-Maguin')
        self.__sg_schoenflies = load_value('str', 'sg_schoenflies', self, valuerequired=True,
                                           modelpath='space-group.Schoenflies')
        self.__crystalfamily = load_value('str', 'crystalfamily', self, valuerequired=True,
                                          modelpath='system-info.cell.crystal-family')
        self.__natypes = load_value('int', 'natypes', self, valuerequired=True,
                                    modelpath='system-info.cell.natypes')

        value_objects.extend([
            self.__key, self.__id, self.__url, self.__commonname, self.__prototype,
            self.__pearson, self.__strukturbericht, self.__sg_number, self.__sg_hm,
            self.__sg_schoenflies, self.__crystalfamily, self.__natypes])

        return value_objects


    def load_model(self,
                   model: Union[str, io.IOBase, DM],
                   name: Optional[str] = None):
        """
        Loads record contents from a given model.

        Parameters
        ----------
        model : str or DataModelDict
            The model contents of the record to load.
        name : str, optional
            The name to assign to the record.  Often inferred from other
            attributes if not given.
        """
        super().load_model(model, name=name)

        self.__ucell = None

        # Set name as id if no name given
        try:
            self.name
        except AttributeError:
            self.name = self.id

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
    def commonname(self) -> str:
        """str : A common name associated with the prototype"""
        return self.__commonname.value
    
    @property
    def prototype(self) -> str:
        """str : A prototype composition associated with the prototype"""
        return self.__prototype.value
    
    @property
    def pearson(self) -> str:
        """str : The prototype's Pearson symbol"""
        return self.__pearson.value

    @property
    def strukturbericht(self) -> str:
        """str : The prototype's Strukturbericht symbol"""
        return self.__strukturbericht.value

    @property
    def sg_number(self) -> int:
        """int : The prototype's space group number"""
        return self.__sg_number.value

    @property
    def sg_hm(self) -> str:
        """str : The prototype's space group international symbol"""
        return self.__sg_hm.value

    @property
    def sg_schoenflies(self) -> str:
        """str : The prototype's space group Schoenflies symbol"""
        return self.__sg_schoenflies.value

    @property
    def crystalfamily(self) -> str:
        """str : The prototype's system family"""
        return self.__crystalfamily.value

    @property
    def natypes(self) -> int:
        """int : Number of atom types"""
        return self.__natypes.value

    @property
    def ucell(self) -> System:
        """atomman.System : The unit cell for the prototype""" 
        if self.model is None:
            raise AttributeError('No model information loaded')
        if self.__ucell is None:
            self.__ucell = System(model=self.model)
        return self.__ucell
