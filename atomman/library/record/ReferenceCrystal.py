# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Union, Tuple
import uuid

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_query, load_value

# atomman imports
from ... import System

class ReferenceCrystal(Record):
    """
    Class for representing reference_crystal records that provide the structure
    information for DFT relaxed crystal structures obtained from DFT databases.
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

        # Call super init
        super().__init__(model=model, name=name, database=database)

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'reference_crystal'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'reference-crystal'

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
        self.__sourcename = load_value('str', 'sourcename', self, valuerequired=True,
                                       modelpath='source.name')
        self.__sourcelink = load_value('str', 'sourcelink', self, valuerequired=True,
                                       modelpath='source.link')
        self.__symbols = load_value('strlist', 'symbols', self, valuerequired=True,
                                    modelpath='system-info.symbol')
        self.__composition = load_value('str', 'composition', self, valuerequired=True,
                                       modelpath='system-info.composition')
        self.__crystalfamily = load_value('str', 'crystalfamily', self,
                                       modelpath='system-info.cell.crystal-family')
        self.__natypes = load_value('int', 'natypes', self, valuerequired=True,
                                    modelpath='system-info.cell.natypes')
        self.__a = load_value('float', 'a', self, valuerequired=True,
                              modelpath='system-info.cell.a')
        self.__b = load_value('float', 'b', self, valuerequired=True,
                              modelpath='system-info.cell.b')
        self.__c = load_value('float', 'c', self, valuerequired=True,
                              modelpath='system-info.cell.c')
        self.__alpha = load_value('float', 'alpha', self, valuerequired=True,
                              modelpath='system-info.cell.alpha')
        self.__beta = load_value('float', 'beta', self, valuerequired=True,
                              modelpath='system-info.cell.beta')
        self.__gamma = load_value('float', 'gamma', self, valuerequired=True,
                              modelpath='system-info.cell.gamma')
        self.__ucell = load_value('system_model', 'ucell', self, valuerequired=True,
                                  modelpath="atomic-system")

        value_objects.extend([
            self.__key, self.__id, self.__url, self.__sourcename, self.__sourcelink,
            self.__symbols, self.__composition, self.__crystalfamily, self.__natypes, 
            self.__a, self.__b, self.__c, self.__alpha, self.__beta, self.__gamma,
            self.__ucell])

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
    def sourcename(self) -> str:
        """str : Name of the crystal's source database"""
        return self.__sourcename.value

    @sourcename.setter
    def sourcename(self, val: str):
        self.__sourcename.value = val

    @property
    def sourcelink(self) -> str:
        """str : URL for the crystal's source database"""
        return self.__sourcelink.value

    @sourcelink.setter
    def sourcelink(self, val: str):
        self.__sourcelink.value = val

    @property
    def symbols(self) -> list:
        """list : The list of element model symbols"""
        return self.__symbols.value
    
    @symbols.setter
    def symbols(self, val: Union[str, list]):
        self.__symbols.value = val

    @property
    def composition(self) -> str:
        """str : The crystal's composition"""
        return self.__composition.value
    
    @composition.setter
    def composition(self, val: str):
        self.__composition.value = val

    @property
    def crystalfamily(self) -> str:
        """str : The crystal's system family"""
        return self.__crystalfamily.value
    
    @crystalfamily.setter
    def crystalfamily(self, val: str):
        self.__crystalfamily.value = val

    @property
    def natypes(self) -> int:
        """int : The number of atom types in the unit cell"""
        return self.__natypes.value
    
    @natypes.setter
    def natypes(self, val: int):
        self.__natypes.value = val

    @property
    def a(self) -> float:
        """float : The unit cell's a lattice parameter"""
        return self.__a.value
    
    @a.setter
    def a(self, val: float):
        self.__a.value = val

    @property
    def b(self) -> float:
        """float : The unit cell's b lattice parameter"""
        return self.__b.value
    
    @b.setter
    def b(self, val: float):
        self.__b.value = val

    @property
    def c(self) -> float:
        """float : The unit cell's c lattice parameter"""
        return self.__c.value
    
    @c.setter
    def c(self, val: float):
        self.__c.value = val

    @property
    def alpha(self) -> float:
        """float : The unit cell's alpha lattice angle"""
        return self.__alpha.value
    
    @alpha.setter
    def alpha(self, val: float):
        self.__alpha.value = val

    @property
    def beta(self) -> float:
        """float : The unit cell's beta lattice angle"""
        return self.__beta.value
    
    @beta.setter
    def beta(self, val: float):
        self.__beta.value = val

    @property
    def gamma(self) -> float:
        """float : The unit cell's gamma lattice angle"""
        return self.__gamma.value
    
    @gamma.setter
    def gamma(self, val: float):
        self.__gamma.value = val

    @property
    def ucell(self) -> System:
        """atomman.System : The unit cell system for the crystal"""
        return self.__ucell.value

    @ucell.setter
    def ucell(self, val: System):
        self.__ucell.value = val
        
    def set_ucell_attributes(self):
        """
        auto sets the symbols, composition, crystalfamily, natypes, a, b, c,
        alpha, beta, and gamma class attributes based on the current ucell.
        """
        self.symbols = self.ucell.symbols
        self.composition = self.ucell.composition
        self.crystalfamily = self.ucell.box.identifyfamily()
        self.natypes = self.ucell.natypes
        self.a = self.ucell.box.a
        self.b = self.ucell.box.b
        self.c = self.ucell.box.c
        self.alpha = self.ucell.box.alpha
        self.beta = self.ucell.box.beta
        self.gamma = self.ucell.box.gamma
