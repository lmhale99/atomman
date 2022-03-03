
from DataModelDict import DataModelDict as DM

from yabadaba import query

# iprPy imports
from . import Record

class StackingFault(Record):
    """
    Class for representing stacking_fault records, which collect the parameters
    necessary for atomman to generate and evaluate a particular generalized
    stacking fault map.
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
        return 'stacking_fault'

    @property
    def xsd_filename(self):
        """tuple: The module path and file name of the record's xsd schema"""
        return ('atomman.library.xsd', f'{self.style}.xsd')
    
    @property
    def modelroot(self):
        """str: The root element of the content"""
        return 'stacking-fault'
    
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
        meta['family'] = self.family
        meta['hkl'] = self.parameters['hkl']
        meta['a1vect_uvw'] = self.parameters['a1vect_uvw']
        meta['a2vect_uvw'] = self.parameters['a2vect_uvw']
        if 'shiftindex' in self.parameters:
            meta['shiftindex'] = self.parameters['shiftindex']
        meta['cutboxvector'] = self.parameters['cutboxvector']
        
        return meta

    def pandasfilter(self, dataframe, name=None, key=None, id=None,
                     family=None, hkl=None, shiftindex=None,
                     cutboxvector=None):
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
        family : str or list
            Parent prototype/reference id(s) to parse by.
        hkl : str or list
            Space delimited fault plane(s) to parse by.
        shiftindex : int or list
            shiftindex value(s) to parse by.
        cutboxvector : str or list
            cutboxvector value(s) to parse by.
        
        Returns
        -------
        pandas.Series, numpy.NDArray
            Boolean map of matching values
        """
        matches = (
            query.str_match.pandas(dataframe, 'name', name)
            &query.str_match.pandas(dataframe, 'key', key)
            &query.str_match.pandas(dataframe, 'id', id)
            &query.str_match.pandas(dataframe, 'family', family)
            &query.str_match.pandas(dataframe, 'hkl', hkl)
            &query.str_match.pandas(dataframe, 'shiftindex', shiftindex)
            &query.str_match.pandas(dataframe, 'cutboxvector', cutboxvector)
        )
        return matches

    def mongoquery(self, name=None, key=None, id=None,
                   family=None, hkl=None, shiftindex=None,
                   cutboxvector=None):
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
        family : str or list
            Parent prototype/reference id(s) to parse by.
        hkl : str or list
            Space delimited fault plane(s) to parse by.
        shiftindex : int or list
            shiftindex value(s) to parse by.
        cutboxvector : str or list
            cutboxvector value(s) to parse by.
        
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
        query.str_match.mongo(mquery, f'{root}.system-family', family)
        query.str_match.mongo(mquery, f'{root}.calculation-parameter.hkl', hkl)
        query.str_match.mongo(mquery, f'{root}.calculation-parameter.shiftindex', shiftindex)
        query.str_match.mongo(mquery, f'{root}.calculation-parameter.cutboxvector', cutboxvector)
        
        return mquery

    def cdcsquery(self, key=None, id=None, family=None, hkl=None, shiftindex=None,
                  cutboxvector=None):
        """
        Builds a CDCS-style query based on kwargs values for the record style.
        
        Parameters
        ----------
        id : str or list
            The record id(s) to parse by.
        key : str or list
            The record key(s) to parse by.
        family : str or list
            Parent prototype/reference id(s) to parse by.
        hkl : str or list
            Space delimited fault plane(s) to parse by.
        shiftindex : int or list
            shiftindex value(s) to parse by.
        cutboxvector : str or list
            cutboxvector value(s) to parse by.
        
        Returns
        -------
        dict
            The CDCS-style query
        """
        mquery = {}
        root = self.modelroot

        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.id', id)
        query.str_match.mongo(mquery, f'{root}.system-family', family)
        query.str_match.mongo(mquery, f'{root}.calculation-parameter.hkl', hkl)
        query.str_match.mongo(mquery, f'{root}.calculation-parameter.shiftindex', shiftindex)
        query.str_match.mongo(mquery, f'{root}.calculation-parameter.cutboxvector', cutboxvector)

        return mquery