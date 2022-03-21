# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Union, Tuple

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_query

# https://pandas.pydata.org/
import pandas as pd

class PointDefect(Record):
    """
    Class for representing point_defect records, which collect the parameters
    necessary for atomman to generate a particular point defect.
    """
    def __init__(self,
                 model: Union[str, io.IOBase, DM, None] = None,
                 name: Optional[str] = None):
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
    def style(self) -> str:
        """str: The record style"""
        return 'point_defect'

    @property
    def xsd_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsd schema"""
        return ('atomman.library.xsd', f'{self.style}.xsd')

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'point-defect'
    
    @property
    def key(self) -> str:
        """str : A UUID4 key assigned to the record"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__key

    @key.setter
    def key(self, value: str):
        self.__key = str(value)
    
    @property
    def id(self) -> str:
        """str : A unique id assigned to the record"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__id

    @id.setter
    def id(self, value: str):
        self.__id = str(value)

    @property
    def family(self) -> str:
        """str : The prototype/reference id the defect is defined for"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__family

    @family.setter
    def family(self, value: str):
        self.__family = str(value)

    @property
    def parameters(self) -> list:
        """list : Defect parameters for atomman structure generator"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__parameters

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
        content = self.model[self.modelroot]

        self.key = content['key']
        self.id = content['id']
        self.family = content['system-family']
        self.__parameters = []
        for cp in content.aslist('calculation-parameter'):
            self.__parameters.append(dict(cp))

    def build_model(self) -> DM:
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
        for cp in self.parameters:
            content.append('calculation-parameter', DM(cp))

        self._set_model(model)
        return model

    def metadata(self) -> dict:
        """
        Generates a dict of simple metadata values associated with the record.
        Useful for quickly comparing records and for building pandas.DataFrames
        for multiple records of the same style.
        """
        meta = {}
        meta['name'] = self.name
        meta['id'] = self.id
        meta['family'] = self.family
        
        meta['ptd_type'] = []
        meta['pos'] = []
        meta['atype'] = []
        meta['db_vect'] = []
        meta['scale'] = []
        for cp in self.parameters:
            meta['ptd_type'].append(cp.get('ptd_type', None))
            meta['pos'].append(cp.get('pos', None))
            meta['atype'].append(cp.get('atype', None))
            meta['db_vect'].append(cp.get('db_vect', None))
            meta['scale'].append(cp.get('scale', None))
        
        return meta

    @property
    def queries(self) -> dict:
        """dict: Query objects and their associated parameter names."""
        return {
            'key': load_query(
                style='str_match',
                name='key', 
                path=f'{self.modelroot}.key'),
            'id': load_query(
                style='str_match',
                name='id',
                path=f'{self.modelroot}.id'),
            'family': load_query(
                style='str_match',
                name='family',
                path=f'{self.modelroot}.system-family'),
            'ptd_type': load_query(
                style='list_contains',
                name='ptd_type',
                path=f'{self.modelroot}.calculation-parameter.ptd_type'),
        }

    def pandasfilter(self,
                     dataframe: pd.DataFrame,
                     name: Union[str, list, None] = None,
                     key: Union[str, list, None] = None,
                     id: Union[str, list, None] = None,
                     family: Union[str, list, None] = None,
                     ptd_type: Union[str, list, None] = None) -> pd.Series:
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
        ptd_type : str or list
            Point defect type(s) to parse by.
        
        Returns
        -------
        pandas.Series
            Boolean map of matching values
        """
        matches = super().pandasfilter(dataframe, name=name, key=key, id=id,
                                       family=family, ptd_type=ptd_type)
        return matches

    def mongoquery(self,
                   name: Union[str, list, None] = None,
                   key: Union[str, list, None] = None,
                   id: Union[str, list, None] = None,
                   family: Union[str, list, None] = None,
                   ptd_type: Union[str, list, None] = None) -> dict:
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
        ptd_type : str or list
            Point defect type(s) to parse by.
        
        Returns
        -------
        dict
            The Mongo-style query
        """   
        mquery = super().mongoquery(name=name, key=key, id=id,
                                    family=family, ptd_type=ptd_type)
        return mquery

    def cdcsquery(self,
                  key: Union[str, list, None] = None,
                  id: Union[str, list, None] = None,
                  family: Union[str, list, None] = None,
                  ptd_type: Union[str, list, None] = None) -> dict:
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
        ptd_type : str or list
            Point defect type(s) to parse by.
        
        Returns
        -------
        dict
            The CDCS-style query
        """
        mquery = super().cdcsquery(key=key, id=id,
                                    family=family, ptd_type=ptd_type)
        return mquery