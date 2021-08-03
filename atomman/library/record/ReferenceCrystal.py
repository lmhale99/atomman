from copy import deepcopy
import uuid

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

from datamodelbase.record import Record
from datamodelbase import query

from ... import System
from ...tools import crystalsystem

modelroot = 'reference-crystal'

class ReferenceCrystal(Record):
    
    @property
    def style(self):
        """str: The record style"""
        return 'reference_crystal'

    @property
    def modelroot(self):
        """str: The root element of the content"""
        return modelroot

    @property
    def xsd_filename(self):
        return ('atomman.library.xsd', f'{self.style}.xsd')

    @property
    def id(self):
        """str : The unique id for the record"""
        return self.__id

    @id.setter
    def id(self, value):
        if value is None:
            self.__id = None
        else:
            self.__id = str(value)

    @property
    def key(self):
        """str : A UUID4 key assigned to the record"""
        return self.__key

    @key.setter
    def key(self, value):
        if value is None:
            self.__key = str(uuid.uuid4())
        else:
            self.__key = str(value)

    @property
    def sourcename(self):
        """str : Name of the crystal's source database"""
        return self.__sourcename

    @sourcename.setter
    def sourcename(self, value):
        if value is None:
            self.__sourcename = None
        else:
            self.__sourcename = str(value)

    @property
    def sourcelink(self):
        """str : URL for the crystal's source database"""
        return self.__sourcelink    

    @sourcelink.setter
    def sourcelink(self, value):
        if value is None:
            self.__sourcelink = None
        else:
            self.__sourcelink = str(value)

    @property
    def ucell(self):
        """atomman.System : The unit cell system for the crystal"""
        if self.__ucell is None:
            raise ValueError('ucell information not set')
        elif not isinstance(self.__ucell, System):
            self.__ucell = System(model=self.__ucell)
        return self.__ucell

    @ucell.setter
    def ucell(self, value):        
        self.__ucell = value

    @property
    def composition(self):
        """str : The crystal's composition"""
        if self.__composition is None:
            self.__composition = self.ucell.composition
        return self.__composition

    @property
    def symbols(self):
        """list : The list of element model symbols"""
        if self.__symbols is None:
            self.__symbols = self.ucell.symbols
        return self.__symbols

    @property
    def natoms(self):
        """int : The number of atoms in the unit cell"""
        if self.__natoms is None:
            self.__natoms = self.ucell.natoms
        return self.__natoms

    @property
    def natypes(self):
        """int : The number of unique elements in the unit cell"""
        if self.__natypes is None:
            self.__natypes = self.ucell.natypes
        return self.__natypes

    @property
    def crystalfamily(self):
        """str : The crystal's system family"""
        if self.__crystalfamily is None:
            self.__crystalfamily = crystalsystem.identifyfamily(self.ucell.box)
        return self.__crystalfamily

    @property
    def a(self):
        """float : The unit cell's a lattice parameter"""
        if self.__a is None:
            self.__a = self.ucell.box.a
        return self.__a

    @property
    def b(self):
        """float : The unit cell's b lattice parameter"""
        if self.__b is None:
            self.__b = self.ucell.box.b
        return self.__b

    @property
    def c(self):
        """float : The unit cell's c lattice parameter"""
        if self.__c is None:
            self.__c = self.ucell.box.c
        return self.__c

    @property
    def alpha(self):
        """float : The unit cell's alpha lattice angle"""
        if self.__alpha is None:
            self.__alpha = self.ucell.box.alpha
        return self.__alpha

    @property
    def beta(self):
        """float : The unit cell's beta lattice angle"""
        if self.__beta is None:
            self.__beta = self.ucell.box.beta
        return self.__beta

    @property
    def gamma(self):
        """float : The unit cell's gamma lattice angle"""
        if self.__gamma is None:
            self.__gamma = self.ucell.box.gamma
        return self.__gamma

    def set_values(self, name=None, id=None, key=None, sourcename=None,
                   sourcelink=None, ucell=None):
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

    def build_model(self):
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
        
    def load_model(self, model, name=None):
        super().load_model(model, name=name)        
        crystal = self.model[modelroot]
        
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
        except:
            self.name = self.id

    def metadata(self):
        """
        Converts the structured content to a simpler dictionary.
        
        Returns
        -------
        dict
            A dictionary representation of the record's content.
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

    @staticmethod
    def pandasfilter(dataframe, name=None, key=None,
                     id=None, sourcename=None,
                     sourcelink=None, crystalfamily=None, composition=None,
                     symbols=None, natoms=None, natypes=None):

        matches = (
            query.str_match.pandas(dataframe, 'name', name)
            &query.str_match.pandas(dataframe, 'key', key)
            &query.str_match.pandas(dataframe, 'id', id)
            &query.str_match.pandas(dataframe, 'sourcename', sourcename)
            &query.str_match.pandas(dataframe, 'sourcelink', sourcelink)
            &query.str_match.pandas(dataframe, 'crystalfamily', crystalfamily)
            &query.str_match.pandas(dataframe, 'composition', composition)
            &query.in_list.pandas(dataframe, 'symbols', symbols)
            &query.str_match.pandas(dataframe, 'natoms', natoms)
            &query.str_match.pandas(dataframe, 'natypes', natypes)
        )
        return matches

    @staticmethod
    def mongoquery(name=None, key=None,
                   id=None, sourcename=None,
                   sourcelink=None, crystalfamily=None, composition=None,
                   symbols=None, natoms=None, natypes=None):

        mquery = {}
        query.str_match.mongo(mquery, f'name', name)
        root = f'content.{modelroot}'

        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.id', id)
        query.str_match.mongo(mquery, f'{root}.source.name', sourcename)
        query.str_match.mongo(mquery, f'{root}.sourcelink', sourcelink)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.crystal-family', crystalfamily)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.composition', composition)
        query.in_list.mongo(mquery, f'{root}.system-info.symbol', symbols)
        query.str_match.mongo(mquery, f'{root}.atomic-system.atoms.natoms', natoms)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.natypes', natypes)

        return mquery

    @staticmethod
    def cdcsquery(key=None, id=None, sourcename=None,
                  sourcelink=None, crystalfamily=None, composition=None,
                  symbols=None, natoms=None, natypes=None):

        mquery = {}
        root = modelroot
        
        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.id', id)
        query.str_match.mongo(mquery, f'{root}.source.name', sourcename)
        query.str_match.mongo(mquery, f'{root}.sourcelink', sourcelink)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.crystal-family', crystalfamily)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.composition', composition)
        query.in_list.mongo(mquery, f'{root}.system-info.symbol', symbols)
        query.str_match.mongo(mquery, f'{root}.atomic-system.atoms.natoms', natoms)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.natypes', natypes)

        return mquery