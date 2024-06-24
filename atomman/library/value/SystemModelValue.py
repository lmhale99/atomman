from typing import Optional, Union, Any

from DataModelDict import DataModelDict as DM

import numpy as np

from ... import System
from ...load import load

from yabadaba.value import Value
from yabadaba import load_query
from yabadaba.record import Record

class SystemModelValue(Value):
    
    def __init__(self,
                 name: str,
                 record: Record,
                 defaultvalue: Optional[Any] = None,
                 valuerequired: bool = False,
                 allowedvalues: Optional[tuple] = None,
                 metadatakey: Union[str, bool, None] = None,
                 metadataparent: Optional[str] = None,
                 modelpath: str = 'system-model'):
        """
        Initialize a general Parameter object.

        Parameters
        ----------
        name : str
            The name of the parameter.  This should correspond to the name of
            the associated class attribute.
        record : Record
            The Record object that the Parameter is used with.
        defaultvalue : any or None, optional
            The default value to use for the property.  The default value of
            None indicates that there is no default value.
        valuerequired: bool
            Indicates if a value must be given for the property.  If True, then
            checks will be performed that a value is assigned to the property.
        allowedvalues : tuple or None, optional
            A list/tuple of values that the parameter is restricted to have.
            Setting this to None (default) indicates any value is allowed.
        metadatakey: str, bool or None, optional
            The key name to use for the property when constructing the record
            metadata dict.  If set to None (default) then name will be used for
            metadatakey.  If set to False then the parameter will not be
            included in the metadata dict.
        metadataparent: str or None, optional
            If given, then this indicates that the metadatakey is actually an
            element of a dict in metadata with this name.  This allows for limited
            support for metadata having embedded dicts.
        modelpath: str, optional
            The period-delimited path after the record root element for
            where the parameter will be found in the built data model.  If set
            to None (default) then 'system-model' is used.
        """

        super().__init__(name, record, defaultvalue=defaultvalue,
                         valuerequired=valuerequired, allowedvalues=allowedvalues,
                         metadatakey=metadatakey, metadataparent=metadataparent,
                         modelpath=modelpath)

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