from copy import deepcopy
from DataModelDict import DataModelDict as DM

from datamodelbase.record import Record
from datamodelbase import query

from ... import System

modelroot = 'crystal-prototype'

class CrystalPrototype(Record):
    
    @property
    def style(self):
        """str: The record style"""
        return 'crystal_prototype'

    @property
    def modelroot(self):
        """str: The root element of the content"""
        return modelroot

    @property
    def xsd_filename(self):
        return ('atomman.library.xsd', 'crystal_prototype.xsd')

    def load_model(self, model, name=None):
        super().load_model(model, name=name)        
        proto = self.model[modelroot]
        
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
        except:
            self.name = self.id

    @property
    def id(self):
        return self.__id

    @property
    def key(self):
        return self.__key

    @property
    def commonname(self):
        return self.__commonname
    
    @property
    def prototype(self):
        return self.__prototype
    
    @property
    def pearson(self):
        return self.__pearson

    @property
    def strukturbericht(self):
        return self.__strukturbericht

    @property
    def sg_number(self):
        return self.__sg_number

    @property
    def sg_hm(self):
        return self.__sg_hm

    @property
    def sg_schoenflies(self):
        return self.__sg_schoenflies

    @property
    def crystalfamily(self):
        return self.__crystalfamily

    @property
    def natypes(self):
        return self.__natypes

    @property
    def ucell(self):
        if self.__ucell is None:
            self.__ucell = System(model=self.model)
        return self.__ucell

    def build_model(self):
        return deepcopy(self.model)

    def metadata(self):
        """Returns the simple metadata fields of the record as a flat dict."""
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

    @staticmethod
    def pandasfilter(dataframe, name=None, id=None,
                     key=None, commonname=None, prototype=None,
                     pearson=None, strukturbericht=None, sg_number=None,
                     sg_hm=None, sg_schoenflies=None, crystalfamily=None,
                     natypes=None):
        matches = (
            query.str_match.pandas(dataframe, 'name', name)
            &query.str_match.pandas(dataframe, 'key', key)
            &query.str_match.pandas(dataframe, 'id', id)
            &query.str_match.pandas(dataframe, 'commonname', commonname)
            &query.str_match.pandas(dataframe, 'prototype', prototype)
            &query.str_match.pandas(dataframe, 'pearson', pearson)
            &query.str_match.pandas(dataframe, 'strukturbericht', strukturbericht)
            &query.str_match.pandas(dataframe, 'sg_number', sg_number)
            &query.str_match.pandas(dataframe, 'sg_hm', sg_hm)
            &query.str_match.pandas(dataframe, 'sg_schoenflies', sg_schoenflies)
            &query.str_match.pandas(dataframe, 'crystalfamily', crystalfamily)
            &query.str_match.pandas(dataframe, 'natypes', natypes)
        )
        return matches

    @staticmethod
    def mongoquery(name=None, id=None, key=None, commonname=None,
                   prototype=None, pearson=None, strukturbericht=None,
                   sg_number=None, sg_hm=None, sg_schoenflies=None,
                   crystalfamily=None, natypes=None):

        mquery = {}
        query.str_match.mongo(mquery, f'name', name)
        root = f'content.{modelroot}'

        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.id', id)
        query.str_match.mongo(mquery, f'{root}.name', commonname)
        query.str_match.mongo(mquery, f'{root}.prototype', prototype)
        query.str_match.mongo(mquery, f'{root}.Pearson-symbol', pearson)
        query.str_match.mongo(mquery, f'{root}.Strukturbericht', strukturbericht)
        query.str_match.mongo(mquery, f'{root}.space-group.number', sg_number)
        query.str_match.mongo(mquery, f'{root}.space-group.Hermann-Maguin', sg_hm)
        query.str_match.mongo(mquery, f'{root}.space-group.Schoenflies', sg_schoenflies)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.crystal-family', crystalfamily)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.natypes', natypes)

        return mquery

    @staticmethod
    def cdcsquery(id=None, key=None, commonname=None,
                  prototype=None, pearson=None, strukturbericht=None,
                  sg_number=None, sg_hm=None, sg_schoenflies=None,
                  crystalfamily=None, natypes=None):

        mquery = {}
        root = modelroot

        query.str_match.mongo(mquery, f'{root}.key', key)
        query.str_match.mongo(mquery, f'{root}.id', id)
        query.str_match.mongo(mquery, f'{root}.name', commonname)
        query.str_match.mongo(mquery, f'{root}.prototype', prototype)
        query.str_match.mongo(mquery, f'{root}.Pearson-symbol', pearson)
        query.str_match.mongo(mquery, f'{root}.Strukturbericht', strukturbericht)
        query.str_match.mongo(mquery, f'{root}.space-group.number', sg_number)
        query.str_match.mongo(mquery, f'{root}.space-group.Hermann-Maguin', sg_hm)
        query.str_match.mongo(mquery, f'{root}.space-group.Schoenflies', sg_schoenflies)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.crystal-family', crystalfamily)
        query.str_match.mongo(mquery, f'{root}.system-info.cell.natypes', natypes)
        
        return mquery