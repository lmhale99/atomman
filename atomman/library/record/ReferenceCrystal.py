# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Union, Tuple
import uuid

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_query

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
        self.__composition = None
        self.__symbols = None
        self.__natoms = None
        self.__natypes = None
        self.__crystalfamily = None
        self.__a = None
        self.__b = None
        self.__c = None
        self.__alpha = None
        self.__beta = None
        self.__gamma = None
        self.__ucell = None

        # Call super init
        super().__init__(model=model, name=name, database=database)

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
    
    @property
    def id(self) -> str:
        """str : The unique id for the record"""
        return self.__id

    @id.setter
    def id(self, value: str):
        if value is None:
            self.__id = None
        else:
            self.__id = str(value)

    @property
    def key(self) -> str:
        """str : A UUID4 key assigned to the record"""
        return self.__key

    @key.setter
    def key(self, value: str):
        if value is None:
            self.__key = str(uuid.uuid4())
        else:
            self.__key = str(value)

    @property
    def sourcename(self) -> str:
        """str : Name of the crystal's source database"""
        return self.__sourcename

    @sourcename.setter
    def sourcename(self, value: str):
        if value is None:
            self.__sourcename = None
        else:
            self.__sourcename = str(value)

    @property
    def sourcelink(self) -> str:
        """str : URL for the crystal's source database"""
        return self.__sourcelink    

    @sourcelink.setter
    def sourcelink(self, value: str):
        if value is None:
            self.__sourcelink = None
        else:
            self.__sourcelink = str(value)

    @property
    def ucell(self) -> System:
        """atomman.System : The unit cell system for the crystal"""
        if self.__ucell is None:
            raise ValueError('ucell information not set')
        elif not isinstance(self.__ucell, System):
            self.__ucell = System(model=self.__ucell)
        return self.__ucell

    @ucell.setter
    def ucell(self, value: System):
        if isinstance(value, System):
            self.__ucell = value
        else:
            raise TypeError('ucell must be an atomman.System')

    @property
    def composition(self) -> str:
        """str : The crystal's composition"""
        if self.__composition is None:
            self.__composition = self.ucell.composition
        return self.__composition

    @property
    def symbols(self) -> list:
        """list : The list of element model symbols"""
        if self.__symbols is None:
            self.__symbols = self.ucell.symbols
        return self.__symbols

    @property
    def natoms(self) -> int:
        """int : The number of atoms in the unit cell"""
        if self.__natoms is None:
            self.__natoms = self.ucell.natoms
        return self.__natoms

    @property
    def natypes(self) -> int:
        """int : The number of atom types in the unit cell"""
        if self.__natypes is None:
            self.__natypes = self.ucell.natypes
        return self.__natypes

    @property
    def crystalfamily(self) -> str:
        """str : The crystal's system family"""
        if self.__crystalfamily is None:
            self.__crystalfamily = self.ucell.box.identifyfamily()
        return self.__crystalfamily

    @property
    def a(self) -> float:
        """float : The unit cell's a lattice parameter"""
        if self.__a is None:
            self.__a = self.ucell.box.a
        return self.__a

    @property
    def b(self) -> float:
        """float : The unit cell's b lattice parameter"""
        if self.__b is None:
            self.__b = self.ucell.box.b
        return self.__b

    @property
    def c(self) -> float:
        """float : The unit cell's c lattice parameter"""
        if self.__c is None:
            self.__c = self.ucell.box.c
        return self.__c

    @property
    def alpha(self) -> float:
        """float : The unit cell's alpha lattice angle"""
        if self.__alpha is None:
            self.__alpha = self.ucell.box.alpha
        return self.__alpha

    @property
    def beta(self) -> float:
        """float : The unit cell's beta lattice angle"""
        if self.__beta is None:
            self.__beta = self.ucell.box.beta
        return self.__beta

    @property
    def gamma(self) -> float:
        """float : The unit cell's gamma lattice angle"""
        if self.__gamma is None:
            self.__gamma = self.ucell.box.gamma
        return self.__gamma

    def set_values(self,
                   name: Optional[str] = None,
                   id: Optional[str] = None,
                   key: Optional[str] = None,
                   sourcename: Optional[str] = None,
                   sourcelink: Optional[str] = None,
                   ucell: Optional[System] = None):
        """
        Sets multiple object values.

        Parameters
        ----------
        name : str, optional
            The name to use for saving the record.  Either name or id should
            be given as they are treated as aliases for this record style.
        id : str, optional
            The unique identifier for the record.  Should be composed of a
            source database tag plus the source database's unique identifier.
        key : str, optional
            A UUID4 key assigned to the record.  Note that if not given a new
            random key will be assigned and therefore the keys might not match
            if similar records were generated independently.
        sourcename : str, optional
            The name of the source database where the reference record was
            retrieved.
        sourcelink : str, optional
            The URL to the source database where the reference record was
            retrieved.
        ucell : atomman.System, optional
            A small unit cell system associated with the reference crystal.
        """
        
        if name is None and id is not None:
            self.name = id
            self.id = id
        elif name is not None and id is None:
            self.name = name
            self.id = name
        else:
            self.name = name
            self.id = id

        self.key = key
        self.sourcename = sourcename
        self.sourcelink = sourcelink
        if ucell is not None:
            self.ucell = ucell
        
        self.__symbols = None
        self.__composition = None
        self.__crystalfamily = None
        self.__natypes = None
        self.__a = None
        self.__b = None
        self.__c = None
        self.__alpha = None
        self.__beta = None
        self.__gamma = None

    def build_model(self) -> DM:
        """
        Returns the object info as data model content
        
        Returns
        ----------
        DataModelDict
            The data model content.
        """

        refmodel = DM()
        
        refmodel['key'] = self.key
        refmodel['id'] = self.id

        refmodel['source'] = DM()
        refmodel['source']['name'] = self.sourcename
        refmodel['source']['link'] = self.sourcelink
        
        refmodel['system-info'] = DM()
        symbols = self.symbols
        if len(symbols) == 1:
            refmodel['system-info']['symbol'] = self.symbols[0]
        else:
            refmodel['system-info']['symbol'] = self.symbols
        refmodel['system-info']['composition'] = self.composition
        
        refmodel['system-info']['cell'] = DM()
        refmodel['system-info']['cell']['crystal-family'] = self.crystalfamily
        refmodel['system-info']['cell']['natypes'] = self.natypes
        refmodel['system-info']['cell']['a'] = self.a
        refmodel['system-info']['cell']['b'] = self.b
        refmodel['system-info']['cell']['c'] = self.c
        refmodel['system-info']['cell']['alpha'] = self.alpha
        refmodel['system-info']['cell']['beta'] = self.beta
        refmodel['system-info']['cell']['gamma'] = self.gamma
        
        refmodel['atomic-system'] = self.ucell.model()['atomic-system']

        model = DM([('reference-crystal', refmodel)])
        self._set_model(model)
        return model
        
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
        crystal = self.model[self.modelroot]
        
        self.__key = crystal['key']
        self.__id = crystal['id']
        self.__sourcename = crystal['source']['name']
        self.__sourcelink = crystal['source']['link']
        
        self.__symbols = crystal['system-info'].aslist('symbol')
        self.__composition = crystal['system-info']['composition']
        self.__crystalfamily = crystal['system-info']['cell']['crystal-family']
        self.__natypes = crystal['system-info']['cell']['natypes']
        self.__a = crystal['system-info']['cell']['a']
        self.__b = crystal['system-info']['cell']['b']
        self.__c = crystal['system-info']['cell']['c']
        self.__alpha = crystal['system-info']['cell']['alpha']
        self.__beta = crystal['system-info']['cell']['beta']
        self.__gamma = crystal['system-info']['cell']['gamma']

        self.__natoms = crystal['atomic-system']['atoms']['natoms']

        self.__ucell = crystal

        # Set name as id if no name given
        try:
            self.name
        except AttributeError:
            self.name = self.id

    def metadata(self) -> dict:
        """
        Generates a dict of simple metadata values associated with the record.
        Useful for quickly comparing records and for building pandas.DataFrames
        for multiple records of the same style.
        """
        params = {}
        params['name'] = self.name
        params['key'] = self.key
        params['id'] = self.id
        params['sourcename'] = self.sourcename
        params['sourcelink'] = self.sourcelink

        params['crystalfamily'] = self.crystalfamily
        params['natypes'] = self.natypes
        params['symbols'] = self.symbols
        params['composition'] = self.composition

        params['a'] = self.a
        params['b'] = self.b
        params['c'] = self.c
        params['alpha'] = self.alpha
        params['beta'] = self.beta
        params['gamma'] = self.gamma
        params['natoms'] = self.natoms

        return params

    @property
    def queries(self) -> dict:
        """dict: Query objects and their associated parameter names."""
        return {
            'key': load_query(
                style='str_match',
                name='key', 
                path=f'{self.modelroot}.key',
                description="search by reference crystal's UUID key"),
            'id': load_query(
                style='str_match',
                name='id',
                path=f'{self.modelroot}.id',
                description="search by reference crystal's id"),
            'sourcename': load_query(
                style='str_match',
                name='sourcename',
                path=f'{self.modelroot}.source.name',
                description='search by name of the source'),
            'sourcelink': load_query(
                style='str_match',
                name='sourcelink',
                path=f'{self.modelroot}.source.link',
                description='search by URL host link for the source'),
            'crystalfamily': load_query(
                style='str_match',
                name='crystalfamily',
                path=f'{self.modelroot}.system-info.cell.crystal-family',
                description="search by reference crystal's crystal family"),
            'composition': load_query(
                style='str_match',
                name='composition',
                path=f'{self.modelroot}.system-info.composition',
                description="search by reference crystal's composition"),
            'symbols': load_query(
                style='list_contains',
                name='symbols',
                path=f'{self.modelroot}.system-info.symbol',
                description="search by reference crystal's symbols"),
            'natoms': load_query(
                style='int_match',
                name='natoms',
                path=f'{self.modelroot}.atomic-system.atoms.natoms',
                description="search by number of atoms in the reference crystal"),
            'natypes': load_query(
                style='int_match',
                name='natypes',
                path=f'{self.modelroot}.system-info.cell.natypes',
                description="search by number of atom types in the reference crystal"),
        }
