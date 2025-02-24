from typing import Optional, Union, Any

import numpy as np

from yabadaba.value import Value
from yabadaba import load_query

class VectorValue(Value):
    
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
                 allowedlengths: Union[int, tuple, None] = None):
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
        allowedlengths: int, tuple, or None, optional
            The length(s) that the vector can have.  The default value of None
            will perform no checks.
        """
        # Handle allowedlengths values
        if allowedlengths is None:
            self.__allowedlengths = None
        elif isinstance(allowedlengths, int):
            self.__allowedlengths = (allowedlengths,)
        else:
            allowedlengths = tuple(allowedlengths)
            for v in allowedlengths:
                if not isinstance(v, int):
                    raise TypeError('allowedlengths must be int(s)')
            self.__allowedlengths = allowedlengths

        super().__init__(name, record, defaultvalue=defaultvalue,
                         valuerequired=valuerequired, allowedvalues=allowedvalues,
                         metadatakey=metadatakey, metadataparent=metadataparent,
                         modelpath=modelpath, description=description)
    
    @property
    def allowedlengths(self):
        """tuple: The lengths the vector value can have"""
        return self.__allowedlengths

    @property
    def str(self) -> Optional[str]:
        """str or None: The Miller string representation of the vector."""
        if self.value is None:
            return self.value
        else:
            return ' '.join([str(i) for i in self.value])
        
    def set_value_mod(self, val):
        if val is None:
            return None
        
        # Parse string
        if isinstance(val, str):
            val = np.asarray(val.split(), dtype=float)
        else:
            val = np.asarray(val, dtype=float)

        if val.ndim != 1:
            raise ValueError(f'{self.name} must be a vector, not multi-dimensional')
        if self.allowedlengths is not None and val.shape[0] not in self.allowedlengths:
            raise ValueError(f'{self.name} is limited to lengths of {self.allowedlengths}')
        return val
    
    def build_model_value(self):
        """Function to modify how values are represented in the model"""
        return self.str
        
    def metadata_value(self):
        """Function to modify how values are represented in the metadata"""
        return self.str