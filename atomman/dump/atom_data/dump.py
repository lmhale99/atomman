# coding: utf-8
# Standard Python libraries
from io import open

# atomman imports
import atomman.unitconvert as uc
from .atoms_prop_info import atoms_prop_info
from .velocities_prop_info import velocities_prop_info
from ...lammps import style
from .. import dump_table

def dump(system, f=None, atom_style=None, units=None, natypes=None,
         potential=None, float_format='%.13f', return_info=True):
    """
    Write a LAMMPS-style atom data file from a System.
    
    Parameters
    ----------
    system : atomman.System 
        The system to write to the atom data file.
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a str.
    atom_style : str, optional
        The LAMMPS atom_style option associated with the data file.  If neither
        atom_style or potential is given, will set atom_style to 'atomic'.
    units : str, optional
        The LAMMPS units option associated with the data file.  If neither
        units or potential is given, will set units 'metal'.
    natypes : int, optional
        Allows the natypes value to be manually changed.  This is needed if
        natypes needs to be greater than the current number of atypes.  If
        neither natypes or potential is given, will use system.natypes.
    potential : atomman.lammps.Potential, optional
        Potential-specific values of atom_style, units, and natypes can be
        extracted from a Potential object.  If both potential and any of the
        individual values are given, the individual values will be used.
    float_format : str, optional
        c-style formatting string for floating point values.  Default value is
        '%.13f'.
    return_info : bool, optional
        Indicates if the LAMMPS command lines associated with reading in the
        file are to be returned as a str.  Default value is True.
    
    Returns
    -------
    content : str
        The data file contents (returned if f is not given).
    read_info : str
        The LAMMPS input command lines to read the created data file in (returned if return_info is True).
    """
    # Wrap atoms because LAMMPS hates atoms out of bounds in atom data files
    system.wrap()

    # Extract potential-based parameters
    if potential is not None:
        if units is None:
            units = potential.units
        if atom_style is None:
            atom_style = potential.atom_style
        if natypes is None:
            natypes = len(potential.normalize_symbols(system.symbols))
    
    # Set default parameter values
    else:
        if units is None:
            units = 'metal'
        if atom_style is None:
            atom_style = 'atomic'
        if natypes is None:
            natypes = system.natypes

    # Get unit information according to the units style
    units_dict = style.unit(units)
    
    # Generate header info
    content = '\n%i atoms\n' % system.natoms
    content += '%i atom types\n' % natypes
    
    # Extract and convert box values
    xlo = uc.get_in_units(system.box.xlo, units_dict['length'])
    xhi = uc.get_in_units(system.box.xhi, units_dict['length'])
    ylo = uc.get_in_units(system.box.ylo, units_dict['length'])
    yhi = uc.get_in_units(system.box.yhi, units_dict['length'])
    zlo = uc.get_in_units(system.box.zlo, units_dict['length'])
    zhi = uc.get_in_units(system.box.zhi, units_dict['length'])
    xy = system.box.xy
    xz = system.box.xz
    yz = system.box.yz
    
    # Write box values
    xf2 = float_format + ' ' + float_format
    content += xf2 % (xlo, xhi) +' xlo xhi\n'
    content += xf2 % (ylo, yhi) +' ylo yhi\n'
    content += xf2 % (zlo, zhi) +' zlo zhi\n'
    if xy != 0.0 or xz != 0.0 or yz != 0.0:
        xf3 = float_format + ' ' + float_format + ' ' + float_format
        content += xf3 % (xy, xz, yz) + ' xy xz yz\n'
    
    # Write atom info
    content += '\nAtoms\n\n'
    prop_info = atoms_prop_info(atom_style, units)
    
    content += dump_table(system, prop_info=prop_info, float_format=float_format)
    
    # Handle velocity information if included
    if 'velocity' in system.atoms_prop():
        
        # Write velocity info
        content += '\nVelocities\n\n'
        prop_info = velocities_prop_info(atom_style, units)
        
        content += dump_table(system, prop_info=prop_info, float_format=float_format)
    
    returns = []
    
    # Save to the file-like object
    if hasattr(f, 'write'):
        f.write(content)
    
    # Save to the file name
    elif f is not None:
        with open(f, 'w') as fp:
            fp.write(content)
    
    # Return as a string
    else:
        returns.append(content)
    
    if return_info is True:
        # Return appropriate units, atom_style, boundary, and read_data LAMMPS commands
        boundary = ''
        for i in range(3):
            if system.pbc[i]:
                boundary += 'p '
            else:
                boundary += 'm '
        
        if isinstance(f, str):
            read_data = 'read_data ' + f
        else:
            read_data = ''
        newline = '\n'
        read_info = newline.join(['# Script and atom data file prepared by atomman package',
                                  '',
                                  'units ' + units,
                                  'atom_style ' + atom_style,
                                  ''
                                  'boundary ' + boundary,
                                  read_data])
        returns.append(read_info)
    
    if len(returns) == 1:
        return returns[0]
    elif len(returns) > 1:
        return tuple(returns)