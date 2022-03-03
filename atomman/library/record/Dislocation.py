from copy import deepcopy

import numpy as np


# iprPy imports
from DataModelDict import DataModelDict as DM

from yabadaba import query

# iprPy imports
from . import Record

class Dislocation(Record):
    """
    Class for representing dislocation records, which collect the parameters
    necessary for atomman to generate a particular dislocation type.
    """
    def __init__(self, model=None, name=None):
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
        """
        if model is not None:
            super().__init__(model=model, name=name)
        elif name is not None:
            self.name = name

    @property
    def style(self):
        """str: The record style"""
        return 'dislocation'

    @property
    def xsd_filename(self):
        """tuple: The module path and file name of the record's xsd schema"""
        return ('atomman.library.xsd', f'{self.style}.xsd')

    @property
    def modelroot(self):
        """str: The root element of the content"""
        return 'dislocation'
    
    @property
    def key(self):
        """str : A UUID4 key assigned to the record"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__key

    @key.setter
    def key(self, value):
        self.__key = str(value)
    
    @property
    def id(self):
        """str : A unique id assigned to the record"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__id

    @id.setter
    def id(self, value):
        self.__id = str(value)

    @property
    def character(self):
        """str : The dislocation's character"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__character

    @character.setter
    def character(self, value):
        self.__character = str(value)

    @property
    def burgers_vector(self):
        """str : String representation of the dislocation's Burgers vector"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__burgers_vector

    @burgers_vector.setter
    def burgers_vector(self, value):
        self.__burgers_vector = str(value)

    @property
    def slip_plane(self):
        """numpy.NDArray : The dislocation's slip plane"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return deepcopy(self.__slip_plane)

    @slip_plane.setter
    def slip_plane(self, value):
        value = np.asarray(value, dtype=int)
        assert value.shape == (3,)
        self.__slip_plane = value

    @property
    def line_direction(self):
        """numpy.NDArray : The dislocation's slip plane"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return deepcopy(self.__line_direction)

    @line_direction.setter
    def line_direction(self, value):
        value = np.asarray(value, dtype=int)
        assert value.shape == (3,)
        self.__line_direction = value

    @property
    def family(self):
        """str : The prototype/reference id the defect is defined for"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__family

    @family.setter
    def family(self, value):
        self.__family = str(value)

    @property
    def parameters(self):
        """dict : Defect parameters for atomman structure generator"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__parameters
    
    def load_model(self, model, name=None):
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
        content = self.model[self.modelroot]

        self.key = content['key']
        self.id = content['id']
        self.character = content['character']
        self.burgers_vector = content['Burgers-vector']
        self.slip_plane = content['slip-plane']
        self.line_direction = content['line-direction']
        self.family = content['system-family']
        self.__parameters = dict(content['calculation-parameter'])

    def build_model(self):
        """
        Returns the object info as data model content
        
        Returns
        ----------
        DataModelDict
            The data model content.
        """
        model = DM()
        model[self.modelroot] = content = DM()

        content['key'] = self.key
        content['id'] = self.id
        content['character'] = self.character
        content['Burgers-vector'] = self.burgers_vector
        content['slip-plane'] = self.slip_plane.tolist()
        content['line-direction'] = self.line_direction.tolist()
        content['system-family'] = self.family
        content['calculation-parameter'] = DM(self.parameters)

        self._set_model(model)
        return model

    def metadata(self):
        """
        Generates a dict of simple metadata values associated with the record.
        Useful for quickly comparing records and for building pandas.DataFrames
        for multiple records of the same style.
        """
        meta = {}
        meta['name'] = self.name
        meta['id'] = self.id
        meta['character'] = self.character
        meta['burgers_vector'] = self.burgers_vector
        meta['slip_plane'] = self.slip_plane
        meta['line_direction'] = self.line_direction
        meta['family'] = self.family
        meta['slip_hkl'] = self.parameters['slip_hkl']
        meta['ξ_uvw'] = self.parameters['ξ_uvw']
        meta['burgers'] = self.parameters['burgers']
        meta['m'] = self.parameters['m']
        meta['n'] = self.parameters['n']
        if 'shift' in self.parameters:
            meta['shift'] = self.parameters['shift']
        if 'shiftscale' in self.parameters:
            meta['shiftscale'] = self.parameters['shiftscale']        
        if 'shiftindex' in self.parameters:
            meta['shiftindex'] = self.parameters['shiftindex']
        
        return meta

    def pandasfilter(self, dataframe, name=None, key=None, id=None,
                     character=None, family=None):
        """
        Filters a pandas.DataFrame based on kwargs values for the record style.
        
        Parameters
        ----------
        dataframe : pandas.DataFrame
            A table of metadata for multiple records of the record style.
        name : str or list
            The record name(s) to parse by.
        id : str or list
            The record id(s) to parse by.
        key : str or list
            The record key(s) to parse by.
        character : str or list
            Dislocation character(s) to parse by.
        family : str or list
            Parent prototype/reference id(s) to parse by.
        
        Returns
        -------
        pandas.Series, numpy.NDArray
            Boolean map of matching values
        """
        matches = (
            query.str_match.pandas(dataframe, 'name', name)
            &query.str_match.pandas(dataframe, 'key', key)
            &query.str_match.pandas(dataframe, 'id', id)
            &query.str_match.pandas(dataframe, 'character', character)
            &query.str_match.pandas(dataframe, 'family', family)
        )
        return matches

    def mongoquery(self, name=None, key=None, id=None, character=None, family=None):
        """
        Builds a Mongo-style query based on kwargs values for the record style.
        
        Parameters
        ----------
        name : str or list
            The record name(s) to parse by.
        id : str or list
            The record id(s) to parse by.
        key : str or list
            The record key(s) to parse by.
        character : str or list
            Dislocation character(s) to parse by.
        family : str or list
            Parent prototype/reference id(s) to parse by.
        
        Returns
        -------
        dict
            The Mongo-style query
        """   
        mquery = {}
        root = f'content.{self.modelroot}'

        query.str_match.mongo(mquery, f'name', name)

        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.id', id)
        query.str_match.mongo(mquery, f'{root}.character', character)
        query.str_match.mongo(mquery, f'{root}.system-family', family)
        
        
        return mquery

    def cdcsquery(self, key=None, id=None, character=None, family=None):
        """
        Builds a CDCS-style query based on kwargs values for the record style.
        
        Parameters
        ----------
        id : str or list
            The record id(s) to parse by.
        key : str or list
            The record key(s) to parse by.
        character : str or list
            Dislocation character(s) to parse by.
        family : str or list
            Parent prototype/reference id(s) to parse by.
        
        Returns
        -------
        dict
            The CDCS-style query
        """
        mquery = {}
        root = self.modelroot

        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.id', id)
        query.str_match.mongo(mquery, f'{root}.character', character)
        query.str_match.mongo(mquery, f'{root}.system-family', family)

        return mquery