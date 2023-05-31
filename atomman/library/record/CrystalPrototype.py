# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Union, Tuple

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_query

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
        proto = self.model[self.modelroot]
        
        self.__key = proto['key']
        self.__id = proto['id']
        self.__commonname = proto['name']
        self.__prototype = proto['prototype']
        self.__pearson = proto['Pearson-symbol']
        self.__strukturbericht = proto['Strukturbericht']
        
        self.__sg_number = proto['space-group']['number']
        self.__sg_hm = proto['space-group']['Hermann-Maguin']
        self.__sg_schoenflies = proto['space-group']['Schoenflies']
        
        self.__crystalfamily = proto['system-info']['cell']['crystal-family']
        self.__natypes = proto['system-info']['cell']['natypes']
        self.__ucell = None

        # Set name as id if no name given
        try:
            self.name
        except AttributeError:
            self.name = self.id

    @property
    def id(self) -> str:
        """str : A unique id assigned to the record"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__id

    @property
    def key(self) -> str:
        """str : A UUID4 key assigned to the record"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__key

    @property
    def commonname(self) -> str:
        """str : A common name associated with the prototype"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__commonname
    
    @property
    def prototype(self) -> str:
        """str : A prototype composition associated with the prototype"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__prototype
    
    @property
    def pearson(self) -> str:
        """str : The prototype's Pearson symbol"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__pearson

    @property
    def strukturbericht(self) -> str:
        """str : The prototype's Strukturbericht symbol"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__strukturbericht

    @property
    def sg_number(self) -> int:
        """int : The prototype's space group number"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__sg_number

    @property
    def sg_hm(self) -> str:
        """str : The prototype's space group international symbol"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__sg_hm

    @property
    def sg_schoenflies(self) -> str:
        """str : The prototype's space group Schoenflies symbol"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__sg_schoenflies

    @property
    def crystalfamily(self) -> str:
        """str : The prototype's system family"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__crystalfamily

    @property
    def natypes(self) -> int:
        """int : Number of atom types"""
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.__natypes

    @property
    def ucell(self) -> System:
        """atomman.System : The unit cell for the prototype""" 
        if self.model is None:
            raise AttributeError('No model information loaded')
        if self.__ucell is None:
            self.__ucell = System(model=self.model)
        return self.__ucell

    def build_model(self) -> DM:
        if self.model is None:
            raise AttributeError('No model information loaded')
        return self.model

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
        params['commonname'] = self.commonname
        params['prototype'] = self.prototype
        params['pearson'] = self.pearson
        params['strukturbericht'] = self.strukturbericht
        params['sg_number'] = self.sg_number
        params['sg_hm'] = self.sg_hm
        params['sg_schoenflies'] = self.sg_schoenflies
        params['crystalfamily'] = self.crystalfamily
        params['natypes'] = self.natypes
        
        return params

    @property
    def queries(self) -> dict:
        """dict: Query objects and their associated parameter names."""
        return {
            'key': load_query(
                style='str_match',
                name='key', 
                path=f'{self.modelroot}.key',
                description="search by crystal prototype's UUID key"),
            'id': load_query(
                style='str_match',
                name='id',
                path=f'{self.modelroot}.id',
                description="search by crystal prototype's id"),
            'commonname': load_query(
                style='str_match',
                name='commonname',
                path=f'{self.modelroot}.name',
                description="search by crystal prototype's common name"),
            'prototype': load_query(
                style='str_match',
                name='prototype',
                path=f'{self.modelroot}.prototype',
                description="search by crystal prototype's prototype composition"),
            'pearson': load_query(
                style='str_match',
                name='pearson',
                path=f'{self.modelroot}.Pearson-symbol',
                description="search by crystal prototype's Pearson symbol"),
            'strukturbericht': load_query(
                style='str_match',
                name='strukturbericht',
                path=f'{self.modelroot}.Strukturbericht',
                description="search by crystal prototype's Strukturbericht symbol"),
            'sg_number': load_query(
                style='int_match',
                name='sg_number',
                path=f'{self.modelroot}.space-group.number',
                description="search by crystal prototype's space group number"),
            'sg_hm': load_query(
                style='str_match',
                name='sg_hm',
                path=f'{self.modelroot}.space-group.Hermann-Maguin',
                description="search by crystal prototype's space group Hermann-Maguin symbol"),
            'sg_schoenflies': load_query(
                style='str_match',
                name='sg_schoenflies',
                path=f'{self.modelroot}.space-group.Schoenflies',
                description="search by crystal prototype's space group Schoenflies symbol"),
            'crystalfamily': load_query(
                style='str_match',
                name='crystalfamily',
                path=f'{self.modelroot}.system-info.cell.crystal-family',
                description="search by crystal prototype's crystal family"),
            'natypes': load_query(
                style='int_match',
                name='natypes',
                path=f'{self.modelroot}.system-info.cell.natypes',
                description="search by number of atom types in the crystal prototype"),
        }
