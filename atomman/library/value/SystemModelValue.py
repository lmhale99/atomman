from typing import Optional, Union, Any

from DataModelDict import DataModelDict as DM

import numpy as np

from ... import System
from ...load import load

from yabadaba.value import Value
from yabadaba import load_query
from yabadaba.record import Record

class SystemModelValue(Value):
    
    @property
    def loadkey(self) -> str:
        """str: The last path name where the content """
        return self.modelpath.split('.')[-1]

    def set_value_mod(self, val: Union[System, DM, None]):
        
        if isinstance(val, System):
            return val
        
        if val is None:
            try:
                val = self.record.model
                if val is None:
                    return None
            except:
                return None
        
        return load('system_model', val, key=self.loadkey)

    def build_model_value(self):
        """Function to modify how values are represented in the model"""
        return self.value.dump('system_model')[self.loadkey]
         
    def load_model_value(self, val):
        return DM([(self.loadkey, val)])
    
    def metadata(self, meta):
        """
        Adds the parameter to the record's metadata dict.

        Parameters
        ----------
        meta : dict
            The metadata dict being built for the record.
        """
        #if self.metadatakey is False:
        #    return
        
        meta['natoms'] = self.value.natoms
        meta['symbols'] = self.value.symbols

    @property
    def _default_queries(self) -> dict:
        """dict: Default query operations to associate with the Parameter style"""
        return {
            'natoms': load_query('int_match',
                                  name='natoms',
                                  parent=self.metadataparent,
                                  path=f'{self.record.modelroot}.{self.modelpath}.atoms.natoms',
                                  description='Return only the records where natoms matches a given value'),
            'symbols': load_query('list_contains',
                                  name='symbols',
                                  path=f'{self.record.modelroot}.{self.modelpath}.atom-type-symbol',
                                  description="Return only the records that contain the listed atomic symbols")
        }