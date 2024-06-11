# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Union, Tuple

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_query, load_value


class PointDefect(Record):
    """
    Class for representing point_defect records, which collect the parameters
    necessary for atomman to generate a particular point defect.
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
        return 'point_defect'

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
        return 'point-defect'
    
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
        
        self.__question = load_value('longstr', 'question', self)
        value_objects.append(self.__question)

        self.__answer = load_value('longstr', 'answer', self)
        value_objects.append(self.__answer)

        return value_objects

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
    def url(self) -> Optional[str]:
        """str : A URL where a copy of the record can be found"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__url

    @url.setter
    def url(self, value: Optional[str]):
        if value is None:
            self.__url = None
        else:
            self.__url = str(value)

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
        self.url = content.get('URL', None)
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
        if self.url is not None:
            content['URL'] = self.url
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
        meta['key'] = self.key
        meta['id'] = self.id
        meta['url'] = self.url
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
                path=f'{self.modelroot}.key',
                description="search by point defect parameter set's UUID key"),
            'id': load_query(
                style='str_match',
                name='id',
                path=f'{self.modelroot}.id',
                description="search by point defect parameter set's id"),
            'family': load_query(
                style='str_match',
                name='family',
                path=f'{self.modelroot}.system-family',
                description="search by the crystal prototype that the point defect parameter set is for"),
            'ptd_type': load_query(
                style='list_contains',
                name='ptd_type',
                path=f'{self.modelroot}.calculation-parameter.ptd_type',
                description="search by the point defect parameter set's defect type"),
        }
