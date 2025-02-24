# coding: utf-8

# Standard Python libraries
import io
from typing import Optional, Union

# http://www.numpy.org/
import numpy as np

# https://pandas.pydata.org/
import pandas as pd

# atomman imports
import atomman.unitconvert as uc
from .process_prop_info import process_prop_info
from ... import Atoms, Box, System
from ...tools import uber_open_rmode

def load(table: Union[str, io.IOBase],
         box: Box,
         symbols: Optional[tuple] = None,
         system: Optional[System] = None,
         prop_name: Optional[list] = None,
         table_name: Optional[list] = None,
         shape: Optional[list] = None,
         unit: Optional[list] = None,
         dtype: Optional[list] = None,
         prop_info: Optional[list] = None,
         skiprows: Optional[int] = None,
         nrows: Optional[int] = None,
         comment: Optional[str] = None,
         usecols: Optional[int] = None,
         header: Union[int, list, str, None] = 'infer') -> System:
    """
    Reads in tabular data into atomic properties.
    
    Parameters
    ----------
    table : str or file-like object
        The table content, file path or file-like object containing the
        content to read.
    box : atomman.Box
        The atomic box to use when generating a System around the data.
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
    system : atomman.System, optional
        The atomic system to load the values to.  If not given, a new system
        will be constructed.
    prop_name : list, optional
         The Atoms properties to generate.  Must be given if prop_info is not.
    table_name : list, optional
        The table column name(s) that correspond to each prop_name.  If not
        given, the table_name values will be based on the prop_name values.
    shape : list, optional
        The shape of each per-atom property.  If not given, will be inferred
        from the length of each table_name value.
    unit : list, optional
        Lists the units for each prop_name as stored in the table.  For a
        value of None, no conversion will be performed for that property.  For
        a value of 'scaled', the corresponding table values will be taken in
        box-scaled units.  If not given, all unit values will be set to None
        (i.e. no conversions).
    dtype : list, optional
        Allows for the data type of each property to be explicitly given.
        Values of None will infer the data type from the corresponding
        property values.  If not given, all values will be None.
    prop_info : list of dict, optional
        Structured form of property conversion parameters, in which each
        dictionary in the list corresponds to a single atoms property.  Each
        dictionary must have a 'prop_name' field, and can optionally have
        'table_name', 'shape', 'unit', and 'dtype' fields.
    skiprows : int, optional
        Number of rows to skip before reading the data.
    usecols : int, optional
        Which columns are to be used. Will be passed to pandas.read_csv() usecols
        option.
    nrows : int, optional
        Number of rows of data to read.
    comment : str, optional
        Specifies a character which indicates all text on a given line after
        is to be considered to be a comment and ignored by parser. This is
        often '#'.
    header : int, list of int, str or None, optional
        Indicates how to handle header info.  Will be passed to pandas.read_csv()
        header option.
        
    Returns
    -------
    atomman.System
        The generated system.
    """
    # Process conversion parameters
    prop_info = process_prop_info(prop_name=prop_name, table_name=table_name,
                                  shape=shape, unit=unit, dtype=dtype, prop_info=prop_info)
    
    # Build list of all table_names
    table_name = []
    for prop in prop_info:
        table_name += prop['table_name']
    
    # Read in table to dataframe
    with uber_open_rmode(table) as f:
        df = pd.read_csv(f, sep=r"\s+", names=table_name, skiprows=skiprows,
                         nrows=nrows, comment=comment, header=header, usecols=usecols)
    if 'id' in df:
        df = df.sort_values('id')
    
    # Generate System
    natoms = len(df)
    if system is None:
        system = System(atoms=Atoms(natoms=natoms), box=box)
    
    # Loop over all properties
    for prop in prop_info:
        pname = prop['prop_name']
        if pname == 'a_id':
            continue
        
        # Get values
        value = df[prop['table_name']].values.reshape((natoms,) + prop['shape'])
        
        if prop['unit'] is not None:
            if prop['unit'] == "scaled":
                value = system.box.position_relative_to_cartesian(value)
            else:
                value = uc.set_in_units(value, prop['unit'])
        value = np.asarray(value, dtype=prop['dtype'])
        system.atoms.view[pname] = value
    
    # Set symbols if given
    if symbols is not None:
        system.symbols = symbols
    
    return system