from copy import deepcopy
from DataModelDict import DataModelDict as DM

from .. import System
from ..tools import identifyfamily
class CrystalPrototype():
    
    def __init__(self, model):
        self.load(model=model)

    def __str__(self):
        return f'Crystal Prototype {self.id}'

    @property
    def id(self):
        return self.__id

    @property
    def key(self):
        return self.__key

    @property
    def name(self):
        return self.__name
    
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
    def ucell(self):
        return self.__ucell

    def load(self, model):
        
        self.__model = model = DM(model)
        
        proto = model['crystal-prototype']
        self.__key = proto['key']
        self.__id = proto['id']
        self.__name = proto['name']
        self.__prototype = proto['prototype']
        self.__pearson = proto['Pearson-symbol']
        self.__strukturbericht = proto['Strukturbericht']
        
        self.__sg_number = proto['space-group']['number']
        self.__sg_hm = proto['space-group']['Hermann-Maguin']
        self.__sg_schoenflies = proto['space-group']['Schoenflies']
        
        self.__ucell = System(model=model)

    def asmodel(self):
        
        return deepcopy(self.__model)

    def asdict(self):
        """
        Converts the structured content to a simpler dictionary.
        
        Returns
        -------
        dict
            A dictionary representation of the record's content.
        """
        params = {}
        params['key'] = self.key
        params['id'] = self.id
        params['name'] = self.name
        params['prototype'] = self.prototype
        params['pearson'] = self.pearson
        params['strukturbericht'] = self.strukturbericht
        
        params['sg_number'] = self.sg_number
        params['sg_hm'] = self.sg_hm
        params['sg_schoenflies'] = self.sg_schoenflies
        
        params['ucell'] = ucell = self.ucell
        params['crystal_family'] = identifyfamily(ucell.box)
        params['natypes'] = ucell.natypes
        
        params['a'] = ucell.box.a
        params['b'] = ucell.box.b
        params['c'] = ucell.box.c
        params['alpha'] = ucell.box.alpha
        params['beta'] = ucell.box.beta
        params['gamma'] = ucell.box.gamma
        params['natoms'] = ucell.natoms
        
        return params