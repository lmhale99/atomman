# coding: utf-8

# Standard Python libraries
import os
import io
from typing import Optional, Union

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

def dump(system,
         f: Union[str, io.IOBase, None] = None,
         box_unit: Optional[list] = None,
         prop_name: Optional[list] = None,
         unit: Optional[list] = None,
         prop_unit: Optional[dict] = None,
         format: Optional[str] = None,
         indent: Optional[int] = None) -> Union[str, DM, None]:
    """
    Dumps a JSON/XML System.model() representation of the system.
    
    Parameters
    ----------
    system : atomman.System
        The system to generate the data model for.
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a DataModelDict.
    box_unit : str, optional
        Length unit to use for the box. Default value is 'angstrom'.
    prop_name : list, optional
        The Atoms properties to include.  If neither prop_name nor prop_unit
        are given, all system properties will be included.
    unit : list, optional
        Lists the units for each prop_name as stored in the table.  For a
        value of None, no conversion will be performed for that property.  For
        a value of 'scaled', the corresponding table values will be taken in
        box-scaled units.  If neither unit nor prop_units given, pos will be
        given in Angstroms and all other values will not be converted.
    prop_unit : dict, optional
        dictionary where the keys are the property keys to include, and
        the values are units to use. If neither unit nor prop_units given,
        pos will be given in Angstroms and all other values will not be
        converted.
    format : str, optional
        File format 'xml' or 'json' to save the content as if f is given.  If
        f is a filename, then the format will be automatically inferred from
        f's extension.  If format is not given and cannot be inferred, then it
        will be set to 'json'.
    indent : int or None, optional
        Indentation option to use for XML/JSON content if f is given.  A value
        of None (default) will add no line separatations or indentations.

    Returns
    -------
    model : DataModelDict.DataModelDict or str
        The generated model representation of the system.  Will be a
        DataModelDict if format is not specified, and a JSON- or XML-formatted
        string if format is specified.  Returned if f is not given.
    """
        
    # Generate model for system
    model = system.model(box_unit=box_unit, prop_name=prop_name, unit=unit,
                         prop_unit=prop_unit)

    # Return DataModelDict or str
    if f is None:
        if format is None:
            return model
        elif format.lower() == 'xml':
            return model.xml(indent=indent)
        elif format.lower() == 'json':
            return model.json(indent=indent)
    
    # Write to file
    else:
        if format is None:
            try:
                format = os.path.splitext(f)[1][1:]
            except:
                format = 'json'
        
        if hasattr(f, 'write'):
            if format.lower() == 'xml':
                return model.xml(fp=f, indent=indent)
            elif format.lower() == 'json':
                return model.json(fp=f, indent=indent)
        
        else:
            with open(f, 'w', encoding='UTF-8') as fp:
                if format.lower() == 'xml':
                    return model.xml(fp=fp, indent=indent)
                elif format.lower() == 'json':
                    return model.json(fp=fp, indent=indent)