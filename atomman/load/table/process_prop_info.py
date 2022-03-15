# coding: utf-8

# Standard Python libraries
from copy import deepcopy
from typing import Optional

# atomman imports
from ...tools import indexstr

def process_prop_info(prop_name: Optional[list] = None,
                      table_name: Optional[list] = None,
                      shape: Optional[list] = None,
                      unit: Optional[list] = None,
                      dtype: Optional[list] = None,
                      prop_info: Optional[list] = None) -> list:
    """
    Handles common setting of prop_info terms
    
    Parameters
    ----------
    prop_name : list, optional
        The Atoms properties to include.  If neither prop_name or prop_info are
        given, all system properties will be included.
    table_name : list, optional
        The dump table column name(s) that correspond to each prop_name.  If not
        given, the table_name values will be based on the prop_name and shape
        values.
    shape : list, optional
        The shape of each per-atom property.  If not given, will be inferred
        from the length of each table_name value.
    unit : list, optional
        Lists the units for each prop_name as stored in the table.  For a
        value of None, no conversion will be performed for that property.  For
        a value of 'scaled', the corresponding table values will be taken in
        box-scaled units.  If not given, unit values will be taken based on
        lammps_units if prop_name corresponds to a standard LAMMPS property,
        otherwise will be set to None (no conversion).
    dtype : list, optional
        Allows for the data type of each property to be explicitly given.
        Values of None will infer the data type from the corresponding
        property values.  If not given, all values will be None.
    prop_info : list of dict, optional
        Structured form of property conversion parameters, in which each
        dictionary in the list corresponds to a single atoms property.  Each
        dictionary must have a 'prop_name' field, and can optionally have
        'table_name', 'shape', 'unit', and 'dtype' fields.
        
    Returns
    -------
    prop_info : list of dict
        The filled-in prop_info structure. Only returned if
        return_prop_info is True.
    """
    
    if prop_info is not None:
        # Check that competing parameters are not given
        try:
            assert prop_name is None
            assert table_name is None
            assert shape is None
            assert unit is None
            assert dtype is None
        except:
            raise ValueError('prop_info cannot be given with prop_name, table_name, shape, unit, or dtype')
        prop_info = deepcopy(prop_info)
    
    # Combine list parameters into prop_info
    else:
        
        # Check that prop_name is given
        if prop_name is None:
            raise ValueError('prop_info or prop_name is required.')
        numprops = len(prop_name)
        
        # Check that sizes match
        try:
            assert table_name is None or len(table_name) == numprops
            assert shape is None or len(shape) == numprops
            assert unit is None or len(unit) == numprops
            assert dtype is None or len(dtype) == numprops
        except:
            raise ValueError('any of prop_name, table_name, shape, unit, and dtype given must be the same length')
        
        # Build prop_info
        prop_info = []
        for i in range(numprops):
            pinfo = {}
            pinfo['prop_name'] = prop_name[i]
            if table_name is not None:
                pinfo['table_name'] = table_name[i]
            if shape is not None:
                pinfo['shape'] = shape[i]
            if unit is not None:
                pinfo['unit'] = unit[i]
            if dtype is not None:
                pinfo['dtype'] = dtype[i]
            prop_info.append(pinfo)
    
    # Fill in missing default values
    for prop in prop_info:
        
        # prop_name is required
        if 'prop_name' not in prop:
            raise KeyError('prop_name required for each entry in prop_info')
        
        # Generate table_name from prop_name and shape
        if 'table_name' not in prop:
            prop['shape'] = prop.get('shape', ())
            prop['table_name'] = []
            for index, istr in indexstr(prop['shape']): # pylint: disable=unused-variable
                prop['table_name'].append(prop['prop_name'] + istr)
        # Make certain table_name is a list
        else:
            try:
                assert not isinstance(prop['table_name'], str)
                prop['table_name'] = list(prop['table_name'])
            except:
                prop['table_name'] = [prop['table_name']]
        
        # Take shape as length of table_name
        if 'shape' not in prop:
            numtnames = len(prop['table_name'])
            if numtnames == 1:
                prop['shape'] = ()
            else:
                prop['shape'] = (numtnames, )
        # Make certain shape is a tuple
        else:
            try:
                prop['shape'] = tuple(prop['shape'])
            except:
                prop['shape'] = (prop['shape'],)
        
        # Default unit and dtype values are None
        prop['unit'] = prop.get('unit', None)
        prop['dtype'] = prop.get('dtype', None)
    
    return prop_info