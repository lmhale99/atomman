from typing import Optional, Union, Any

import numpy as np

from yabadaba.value import Value
from yabadaba import load_query

from ...tools import miller

class MillerValue(Value):
    
    def __init__(self,
                 name: str,
                 record,
                 defaultvalue: Optional[Any] = None,
                 valuerequired: bool = False,
                 allowedvalues: Optional[tuple] = None,
                 metadatakey: Union[str, bool, None] = None,
                 metadataparent: Optional[str] = None,
                 modelpath: Optional[str] = None,
                 description: Optional[str] = None,
                 bracket: str = '[]'):
        """
        Initialize a Value object.

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
            to None (default) then name will be used for modelpath.
        description: str or None, optional
            A short description for the value.  If not given, then the record name
            will be used.
        bracket: str, optional
            The style of Miller vector or plane brackets to use: '[]', '<>', '()',
            or '{}'.  Default value is '[]'.
        """
        self.__bracket = bracket

        super().__init__(name, record, defaultvalue=defaultvalue,
                         valuerequired=valuerequired, allowedvalues=allowedvalues,
                         metadatakey=metadatakey, metadataparent=metadataparent,
                         modelpath=modelpath, description=description)
    
    @property
    def bracket(self) -> str:
        """str: The bracket style to use for the string representation."""
        return self.__bracket
    
    @bracket.setter
    def bracket(self, val:str):
        if val not in ['[]', '<>', '()', '{}']:
            raise ValueError("Invalid bracket style: allowed values are '[]', '<>', '()', and '{}'")

    @property
    def str(self) -> Optional[str]:
        """str or None: The Miller string representation of the vector."""
        if self.value is None:
            return self.value
        else:
            return miller.tostring(self.value, self.bracket)
        
    def set_value_mod(self, val):
        if val is None:
            return None
        
        # Convert val to Miller string for verification
        if not isinstance(val, str):
            val = miller.tostring(val, self.bracket)

        # Convert back to number
        return miller.fromstring(val)
        
    @property
    def _default_queries(self) -> dict:
        """dict: Default query operations to associate with the Parameter style"""
        return {
            self.name: load_query('str_match',
                                  name=self.metadatakey,
                                  parent=self.metadataparent,
                                  path=f'{self.record.modelroot}.{self.modelpath}',
                                  description=f'Return only the records where {self.name} matches a given value')
        }
    
    def build_model_value(self):
        """Function to modify how values are represented in the model"""
        return self.str
        
    def load_model_value(self, val):
        """Function to modify how values are interpreted from the model"""
        return miller.fromstring(val)
        
    def metadata_value(self):
        """Function to modify how values are represented in the metadata"""
        return self.str