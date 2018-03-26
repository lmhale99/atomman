# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)
from io import StringIO, open

# atomman imports
import atomman.unitconvert as uc
from ...lammps import style
from .. import table
from .standard_prop_info import standard_prop_info

def dump(system, f=None, lammps_units='metal', scale=False, prop_name=None,
         table_name=None, shape=None, unit=None, dtype=None,
         prop_info=None, float_format ='%.13f', return_prop_info=False):
    """
    Write a LAMMPS-style atom data file from a System.
    
    Parameters
    ----------
    system : atomman.System 
        The system to write to the atom data file.
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a str.
    lammps_units : str, optional
        The LAMMPS units option associated with the table values.  This is used
        for the box dimensions and default units for standard dump properties 
        (not compute/fix definitions). 
    scale : bool, optional
        Flag indicating if atom positions are to be scaled relative to the box
        (True) or given in absolute Cartesian values (False, default).
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
    float_format : str, optional
        c-style formatting string for floating point values.  Default value is
        '%.13f'.
    return_prop_info : bool, optional
        Flag indicating if the filled-in prop_info is to be returned.  Having
        this allows for 1:1 load/dump conversions.  Default value is False
        (prop_info is not returned).
        
    Returns
    -------
    prop_info : list of dict
        The filled-in prop_info structure (if return_prop_info is True).
    """
    lammps_unit = style.unit(lammps_units)
    all_prop_info = standard_prop_info(lammps_units)
    
    # Generate prop_info if needed
    if prop_info is None:
        
        # Take all atoms properties if no prop_name or prop_info
        if prop_name is None:
            try:
                assert table_name is None 
                assert shape is None 
                assert unit is None
                assert dtype is None
            except:
                raise ValueError('table_name, shape, unit and dtype cannot be given without prop_name')
            prop_name = ['a_id'] + system.atoms_prop()
        
        # Build prop_info
        prop_info = []
        for i in range(len(prop_name)):
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
        
    # Fill in prop_info values
    for prop in prop_info:
        
        if prop['prop_name'] == 'pos' and scale is True:
            match = True
            sprop = {}
            sprop['table_name'] = ['xs', 'ys', 'zs']
            sprop['unit'] = 'scaled'
        else:
            match = False
            for sprop in all_prop_info:
                if prop['prop_name'] == sprop['prop_name']:
                    match = True
                    break
        if match is True:
            for key in ('table_name', 'shape', 'unit', 'dtype'):
                if key in sprop and key not in prop:
                    prop[key] = sprop[key]
        else:
            # Get shapes of properties without table_names
            if 'table_name' not in prop and 'shape' not in prop:
                prop['shape'] = system.atoms.view[prop['prop_name']].shape[1:]
                
    # Process conversion parameters
    prop_info = table.process_prop_info(prop_info=prop_info)
    
    # Open file
    if f is None:
        fp = StringIO()
        fclose = True
        fstr = True
    elif hasattr(f, 'write'):
        fp = f
        fclose = False
        fstr = False
    else:
        fp = open(f, 'w')
        fclose = True
        fstr = False
    
    # Write timestep info
    fp.write('ITEM: TIMESTEP\n')
    try:
        fp.write('%i\n' % system.timestep)
    except:
        fp.write('0\n')
    
    # Write number of atoms
    fp.write('ITEM: NUMBER OF ATOMS\n')
    fp.write('%i\n' % (system.natoms))
    
    # Extract and convert box values
    xlo = uc.get_in_units(system.box.xlo, lammps_unit['length'])
    xhi = uc.get_in_units(system.box.xhi, lammps_unit['length'])
    ylo = uc.get_in_units(system.box.ylo, lammps_unit['length'])
    yhi = uc.get_in_units(system.box.yhi, lammps_unit['length'])
    zlo = uc.get_in_units(system.box.zlo, lammps_unit['length'])
    zhi = uc.get_in_units(system.box.zhi, lammps_unit['length'])
    xy = uc.get_in_units(system.box.xy, lammps_unit['length'])
    xz = uc.get_in_units(system.box.xz, lammps_unit['length'])
    yz = uc.get_in_units(system.box.yz, lammps_unit['length'])
    
    # Compute absolute box bounds
    xlo_bound = xlo + min((0.0, xy, xz, xy + xz))
    xhi_bound = xhi + max((0.0, xy, xz, xy + xz))
    ylo_bound = ylo + min((0.0, yz))
    yhi_bound = yhi + max((0.0, yz))
    zlo_bound = zlo
    zhi_bound = zhi
    
    is_orthogonal = (xy == 0.0 and xz == 0.0 and yz == 0.0)
    
    # Write system boundary info
    if is_orthogonal:
        fp.write('ITEM: BOX BOUNDS')
    else:
        fp.write('ITEM: BOX BOUNDS xy xz yz')
        
    # Write pbc info
    for i in range(3):
        if system.pbc[i]:
            fp.write(' pp')
        else:
            fp.write(' fm')
    fp.write('\n')
        
    # Write system boundary info
    if is_orthogonal:
        xf2 = float_format + ' ' + float_format + '\n'
        fp.write(xf2 % (xlo_bound, xhi_bound))
        fp.write(xf2 % (ylo_bound, yhi_bound))
        fp.write(xf2 % (zlo_bound, zhi_bound))
    else:
        xf3 = float_format + ' ' + float_format + ' ' + float_format + '\n'
        fp.write(xf3 % (xlo_bound, xhi_bound, xy))
        fp.write(xf3 % (ylo_bound, yhi_bound, xz))
        fp.write(xf3 % (zlo_bound, zhi_bound, yz))
    
    # Write atomic header info and prepare outarray for writing
    header = 'ITEM: ATOMS'
    for prop in prop_info:
        header += ' ' + ' '.join(prop['table_name'])
    header += '\n'
    fp.write(header)
    
    table.dump(system, f=fp, prop_info=prop_info, float_format=float_format)
    
    returns = []
    
    if fclose is True:
        if fstr is True:
            returns.append(fp.getvalue())
        fp.close()
    
    if return_prop_info is True:
        returns.append(prop_info)
        
    if len(returns) == 1:
        return returns[0]
    elif len(returns) > 1:
        return tuple(returns)