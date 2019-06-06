# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
import os

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# atomman imports
import atomman.unitconvert as uc
from ...compatibility import range, iteritems, stringtype
from ...tools import identifyfamily

def dump(system, f=None, box_unit=None, prop_name=None, unit=None,
         prop_unit=None, format=None, indent=None, symbols=None):
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
    symbols : list, optional
        list of atom-model symbols corresponding to the atom types.  If not
        given, will use system.symbols.

    Returns
    -------
    model : DataModelDict.DataModelDict or str
        A JSON/XML data model for the current System object.  Returned if f is not given.
    """
        
    # Generate model for system
    model = system.model(box_unit=box_unit, symbols=symbols, prop_name=prop_name,
                         unit=unit, prop_unit=prop_unit)

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
            with open(f, 'w') as fp:
                if format.lower() == 'xml':
                    return model.xml(fp=fp, indent=indent)
                elif format.lower() == 'json':
                    return model.json(fp=fp, indent=indent)