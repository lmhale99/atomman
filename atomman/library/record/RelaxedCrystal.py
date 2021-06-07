from copy import deepcopy
from DataModelDict import DataModelDict as DM

from datamodelbase.record import Record
from datamodelbase import query

import atomman.unitconvert as uc
from ... import System

modelroot = 'relaxed-crystal'

class RelaxedCrystal(Record):
    
    @property
    def style(self):
        """str: The record style"""
        return 'relaxed_crystal'

    @property
    def modelroot(self):
        """str: The root element of the content"""
        return modelroot

    @property
    def xsd_filename(self):
        return ('atomman.library.xsd', 'relaxed_crystal.xsd')

    def load_model(self, model, name=None):
        super().load_model(model, name=name)        
        crystal = self.model[modelroot]
        
        self.__key = crystal['key']
        self.__method = crystal['method']
        self.__standing = crystal['standing']

        self.__potential_LAMMPS_id = crystal['potential-LAMMPS']['id']
        self.__potential_LAMMPS_key = crystal['potential-LAMMPS']['key']
        self.__potential_id = crystal['potential-LAMMPS']['potential']['id']
        self.__potential_key = crystal['potential-LAMMPS']['potential']['key']

        self.__family = crystal['system-info']['family']
        self.__parent_key = crystal['system-info']['parent_key']
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

        self.__potential_energy = uc.value_unit(crystal['potential-energy'])
        self.__cohesive_energy = uc.value_unit(crystal['cohesive-energy'])
        
        self.__ucell = None

        # Set name as key if no name given
        try:
            self.name
        except:
            self.name = self.key

    @property
    def key(self):
        return self.__key

    @property
    def method(self):
        return self.__method
    
    @property
    def standing(self):
        return self.__standing
    
    @property
    def family(self):
        return self.__family

    @property
    def parent_key(self):
        return self.__parent_key

    @property
    def potential_LAMMPS_id(self):
        return self.__potential_LAMMPS_id

    @property
    def potential_LAMMPS_key(self):
        return self.__potential_LAMMPS_key

    @property
    def potential_id(self):
        return self.__potential_id

    @property
    def potential_key(self):
        return self.__potential_key

    @property
    def cohesive_energy(self):
        return self.__cohesive_energy

    @property
    def composition(self):
        return self.__composition

    @property
    def symbols(self):
        return self.__symbols

    @property
    def natoms(self):
        return self.__natoms

    @property
    def natypes(self):
        return self.__natypes

    @property
    def crystalfamily(self):
        return self.__crystalfamily

    @property
    def a(self):
        return self.__a

    @property
    def b(self):
        return self.__b

    @property
    def c(self):
        return self.__c

    @property
    def alpha(self):
        return self.__alpha

    @property
    def beta(self):
        return self.__beta

    @property
    def gamma(self):
        return self.__gamma

    @property
    def ucell(self):
        if self.__ucell is None:
            self.__ucell = System(model=self.model)
        return self.__ucell

    def build_model(self):
        
        return deepcopy(self.model)

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
        params['method'] = self.method
        params['standing'] = self.standing
        params['family'] = self.family
        params['parent_key'] = self.parent_key
        
        params['potential_LAMMPS_id'] = self.potential_LAMMPS_id
        params['potential_LAMMPS_key'] = self.potential_LAMMPS_key
        params['potential_id'] = self.potential_id
        params['potential_key'] = self.potential_key
        
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

        params['potential_energy'] = self.potential_energy
        params['cohesive_energy'] = self.cohesive_energy

        return params

    @staticmethod
    def pandasfilter(dataframe, name=None, key=None,
                     method='dynamic', standing='good',
                     family=None, parent_key=None, 
                     potential_LAMMPS_id=None, potential_LAMMPS_key=None,
                     potential_id=None, potential_key=None,
                     symbols=None, natoms=None, natypes=None):

        matches = (
            query.str_match.pandas(dataframe, 'name', name)
            &query.str_match.pandas(dataframe, 'key', key)
            &query.str_match.pandas(dataframe, 'method', method)
            &query.str_match.pandas(dataframe, 'standing', standing)
            &query.str_match.pandas(dataframe, 'family', family)
            &query.str_match.pandas(dataframe, 'parent_key', parent_key)
            &query.str_match.pandas(dataframe, 'potential_LAMMPS_id', potential_LAMMPS_id)
            &query.str_match.pandas(dataframe, 'potential_LAMMPS_key', potential_LAMMPS_key)
            &query.str_match.pandas(dataframe, 'potential_id', potential_id)
            &query.str_match.pandas(dataframe, 'potential_key', potential_key)
            &query.str_match.pandas(dataframe, 'potential_key', potential_key)
            &query.in_list.pandas(dataframe, 'symbols', symbols)
            &query.str_match.pandas(dataframe, 'natoms', natoms)
            &query.str_match.pandas(dataframe, 'natypes', natypes)
        )
        return matches

    @staticmethod
    def mongoquery(name=None, key=None,
                   method='dynamic', standing='good',
                   family=None, parent_key=None, 
                   potential_LAMMPS_id=None, potential_LAMMPS_key=None,
                   potential_id=None, potential_key=None,
                   symbols=None, natoms=None, natypes=None):

        mquery = {}
        query.str_match.mongo(mquery, f'name', name)
        root = f'content.{modelroot}'

        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.method', method)
        query.str_match.mongo(mquery, f'{root}.standing', standing)
        query.str_match.mongo(mquery, f'{root}.system-info.family', family)
        query.str_match.mongo(mquery, f'{root}.system-info.parent_key', parent_key)
        query.str_match.mongo(mquery, f'{root}.potential-LAMMPS.id', potential_LAMMPS_id)
        query.str_match.mongo(mquery, f'{root}.potential-LAMMPS.key', potential_LAMMPS_key)
        query.str_match.mongo(mquery, f'{root}.potential-LAMMPS.potential.id', potential_id)
        query.str_match.mongo(mquery, f'{root}.potential-LAMMPS.potential.key', potential_key)
        query.in_list.mongo(mquery, f'{root}.system-info.symbol', symbols)
        query.str_match.mongo(mquery, f'{root}.atomic-system.atoms.natoms', natoms)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.natypes', natypes)

        return mquery

    @staticmethod
    def cdcsquery(key=None, method='dynamic', standing='good',
                  family=None, parent_key=None, 
                  potential_LAMMPS_id=None, potential_LAMMPS_key=None,
                  potential_id=None, potential_key=None,
                  symbols=None, natoms=None, natypes=None):

        mquery = {}
        root = modelroot
        
        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.method', method)
        query.str_match.mongo(mquery, f'{root}.standing', standing)
        query.str_match.mongo(mquery, f'{root}.system-info.family', family)
        query.str_match.mongo(mquery, f'{root}.system-info.parent_key', parent_key)
        query.str_match.mongo(mquery, f'{root}.potential-LAMMPS.id', potential_LAMMPS_id)
        query.str_match.mongo(mquery, f'{root}.potential-LAMMPS.key', potential_LAMMPS_key)
        query.str_match.mongo(mquery, f'{root}.potential-LAMMPS.potential.id', potential_id)
        query.str_match.mongo(mquery, f'{root}.potential-LAMMPS.potential.key', potential_key)
        query.in_list.mongo(mquery, f'{root}.system-info.symbol', symbols)
        query.str_match.mongo(mquery, f'{root}.atomic-system.atoms.natoms', natoms)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.natypes', natypes)

        return mquery