# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Tuple, Union

# https://pandas.pydata.org/
import pandas as pd

# atomman imports
import atomman.unitconvert as uc
from .atoms_prop_info import atoms_prop_info
from .velocities_prop_info import velocities_prop_info
from ... import Atoms, Box, System
from ...lammps import style
from .. import load as amload 
from .. import FileFormatError
from ...tools import uber_open_rmode

def load(data: Union[str, io.IOBase],
         pbc: Tuple[bool, bool, bool] = (True, True, True),
         symbols: Optional[tuple] = None,
         atom_style: Optional[str] = None,
         units: str = 'metal') -> System:
    """
    Read a LAMMPS-style atom data file.
    
    Parameters
    ----------
    data : str or file-like object
        The atom data content to read.  Can be str content, path name, or open
        file-like object.
    pbc : list of bool
        Three boolean values indicating which System directions are periodic.
        Default value is (True, True, True).
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
    atom_style : str, optional
        The LAMMPS atom_style option associated with the data file.  Optional
        as the data can list this value in a comment in the Atoms section
        header.  If not given and not found in data, the default value of
        'atomic' is used.
    units : str, optional
        The LAMMPS units option associated with the data file.  Default value
        is 'metal'.
    
    Returns
    -------
    atomman.System
        The corresponding system.  Note all property values will be
        automatically converted to atomman.unitconvert's working units.

    Raises
    ------
    FileNotFoundError
        If data is (likely) a file name and no matching file is found.
    FileFormatError
        If required content fields not found.
    ValueError
        If atom_style is both given as a parameter and found in data, but are
        not the same
    """

    # First pass over file to generate system and locate content
    system, params = firstpass(data, pbc, symbols, units)

    # Set atom_style
    if atom_style is None:
        
        # Use default if no value given/found
        if params['atom_style'] is None:
            atom_style = 'atomic'
        
        # Use value listed in data
        else:
            atom_style = params['atom_style']
    
    # Check that atom_style matches the one in data
    elif params['atom_style'] is not None and atom_style != params['atom_style']:
        raise ValueError(f'given atom_style of {atom_style} differs from value of {params["atom_style"]} found in data')
    
    # Read in Atoms info
    system = read_atoms(data, system, atom_style, units, params['atomsstart'], params['atomscolumns'])
    
    # Read in Velocities info
    system = read_velocities(data, system, atom_style, units, params['velocitiesstart'])
    
    return system

def firstpass(data, pbc, symbols, units):
    """
    Reads through data to extract natoms and box dimensions to construct a
    System sans atomic properties.

    Parameters
    ----------
    data : str or file-like object
        The atom data content to read.  Can be str content, path name, or open
        file-like object.
    pbc : list of bool
        Three boolean values indicating which System directions are periodic.
        Default value is (True, True, True).
    symbols : tuple, optional
        Allows the list of element symbols to be assigned during loading.
    units : str 
        The LAMMPS units option associated with the data file. Default value
        is 'metal'.
    
    Returns
    -------
    system : atomman.System
        The system with the correct number of atoms and no assigned atom
        properties.
    params : dict
        Line indices and flags indicating how to read in the atom properties
        from the file on subsequent passes.

    Raises
    ------
    FileFormatError
        If '# atoms' or box boundaries not found.
    """
    # Get units information
    units_dict = style.unit(units)
    
    # Initialize parameter values
    atomsstart = None
    velocitiesstart = None
    natoms = None
    natypes = None
    firstatoms = False
    atomscolumns = 0
    masses = None
    num_masses_to_read = 0
    xlo = xhi = ylo = yhi = zlo = zhi = None
    xy = 0.0
    xz = 0.0
    yz = 0.0
    i = 0
    
    # Read str and files in the same way
    with uber_open_rmode(data) as fp:
        
        # Loop over all lines in fp
        for i, fullline in enumerate(fp):
            try:
                fullline = fullline.decode('UTF-8')
            except:
                pass
            
            # Remove comments after '#'
            try:
                comment_index = fullline.index('#')
            except:
                line = fullline
            else:
                line = fullline[:comment_index]
            
            terms = line.split()

            # Skip blank lines
            if len(terms)>0:
                
                # Read number of atoms 
                if len(terms) == 2 and terms[1] == 'atoms':
                    natoms = int(terms[0])

                # Read number of atom types
                elif len(terms) == 3 and terms[1] == 'atom' and terms[2] == 'types': 
                    natypes = int(terms[0])
                
                # Read boundary info
                elif len(terms) == 4 and terms[2] == 'xlo' and terms[3] == 'xhi':
                    xlo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    xhi = uc.set_in_units(float(terms[1]), units_dict['length'])

                elif len(terms) == 4 and terms[2] == 'ylo' and terms[3] == 'yhi':
                    ylo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    yhi = uc.set_in_units(float(terms[1]), units_dict['length'])

                elif len(terms) == 4 and terms[2] == 'zlo' and terms[3] == 'zhi':
                    zlo = uc.set_in_units(float(terms[0]), units_dict['length'])
                    zhi = uc.set_in_units(float(terms[1]), units_dict['length'])

                elif len(terms) == 6 and terms[3] == 'xy' and terms[4] == 'xz' and terms[5] == 'yz':
                    xy = uc.set_in_units(float(terms[0]), units_dict['length'])
                    xz = uc.set_in_units(float(terms[1]), units_dict['length'])
                    yz = uc.set_in_units(float(terms[2]), units_dict['length'])
                
                # Identify starting line number for Atoms data
                elif len(terms) == 1 and terms[0] == 'Atoms':
                    atomsstart = i + 1
                    firstatoms = True

                    # Check for atom_style comment
                    try: 
                        comment_index = fullline.index('#')
                    except:
                        atom_style = None
                    else:
                        atom_style = fullline[comment_index + 1:].strip()
                
                # Count number of columns in Atoms table
                elif firstatoms:
                    atomscolumns = len(terms)
                    firstatoms = False
                
                # Identify starting line for Masses data
                elif len(terms) == 1 and terms[0] == 'Masses':
                    if natypes is None:
                        raise FileFormatError('# atom types must appear before Masses list')
                    masses = [None for i in range(natypes)]
                    num_masses_to_read = natypes
                
                # Read masses
                elif num_masses_to_read > 0:
                    read_mass(terms, masses)
                    num_masses_to_read -= 1

                # Identify starting line number for Velocity data
                elif len(terms) == 1 and terms[0] == 'Velocities':
                    velocitiesstart = i + 1
    
    if i == 0:
        raise FileNotFoundError(f'File {data} not found')

    if natoms is None:
        raise FileFormatError('# atoms not found')

    if xlo is None or xhi is None:
        raise FileFormatError('xlo, xhi box dimensions missing')

    if ylo is None or yhi is None:
        raise FileFormatError('ylo, yhi box dimensions missing')

    if zlo is None or zhi is None:
        raise FileFormatError('zlo, zhi box dimensions missing')

    if atomsstart is None:
        raise FileFormatError('Atoms section missing')

    # Create system with natoms
    box = Box(xlo=xlo, xhi=xhi,
              ylo=ylo, yhi=yhi,
              zlo=zlo, zhi=zhi,
              xy=xy, xz=xz, yz=yz)
    atoms = Atoms(natoms=natoms)
    system = System(box=box, atoms=atoms, pbc=pbc, symbols=symbols,
                    masses=masses)

    # Compile dict of params
    params = {}
    params['atomsstart'] = atomsstart
    params['velocitiesstart'] = velocitiesstart
    params['atomscolumns'] = atomscolumns
    params['atom_style'] = atom_style

    return system, params

def read_atoms(data, system, atom_style, units, atomsstart, atomscolumns):
    """
    Reads in an "Atoms" table from data.

    Parameters
    ----------
    data : str or file-like object
        The atom data content to read.  Can be str content, path name, or open
        file-like object.
    system : atomman.System
        The system to add Atoms table info to.
    atom_style :str
        The LAMMPS atom_style option associated with the data file.  Default
        value is 'atomic'.
    units : str 
        The LAMMPS units option associated with the data file. Default value
        is 'metal'.
    atomsstart : int or None
        The line of the file where the Atoms table content starts.
    atomscolumns : int
        How many columns are in the Atoms table.  Used to determine if wrap
        flags are included.

    Returns
    -------
    atomman.System
        The system with atom properties from the Atoms table assigned.
    """
    
    if atomsstart is not None:
        prop_info = atoms_prop_info(atom_style, units)
        ncols = countreadcolumns(prop_info)
        
        # Read Atoms table
        system = amload('table', data, box=system.box, system=system, 
                        prop_info=prop_info, skiprows=atomsstart,
                        nrows=system.natoms, comment='#',
                        header=None, usecols=range(ncols))
        
        # Check if image flags are included
        if atomscolumns == ncols + 3:
            
            # Read image flags
            with uber_open_rmode(data) as f:
                imageflags = pd.read_csv(f, sep=r"\s+", names=['bx', 'by', 'bz'],
                                        skiprows=atomsstart, nrows=system.natoms, comment='#',
                                        header=None, usecols=range(ncols, atomscolumns),
                                        dtype='int64')

            # Wrap atoms to correct images
            shift = imageflags.values.dot(system.box.vects)
            system.atoms.pos[:] += shift
        
        # Check for correct number of columns
        elif ncols != atomscolumns:
            raise FileFormatError(f'atom_style={atom_style} requires {ncols} or {ncols+3} Atoms table columns but {atomscolumns} found')

    return system

def countreadcolumns(prop_info):
        """
        Counts how many columns are expected from prop info
        """
        count = 0
        for prop in prop_info:
            if isinstance(prop['table_name'], str):
                count += 1
            else:
                count += len(prop['table_name'])
        return count

def read_mass(terms, masses):
    """
    Reads in a line in the Masses table.

    Parameters
    ----------
    terms : list
        The space-delimited terms in an input line
    masses : list
        The list of already assigned masses
    """
    # Check that Masses line is correctly formatted
    try:
        assert len(terms) == 2
        atype = int(terms[0])
        assert atype > 0 and atype <= len(masses)
        mass = float(terms[1])
        assert mass > 0
    except:
        raise FileFormatError('Invalid mass term')
    
    if masses[atype - 1] is None:
        masses[atype - 1] = mass
    else:
        raise FileFormatError(f'Multiple masses listed for atom type {atype}')

def read_velocities(data, system, atom_style, units, velocitiesstart):
    """
    Reads in an "Velocities" table from data.

    Parameters
    ----------
    data : str or file-like object
        The atom data content to read.  Can be str content, path name, or open
        file-like object.
    system : atomman.System
        The system to add Velocities table info to.
    atom_style :str
        The LAMMPS atom_style option associated with the data file.  Default
        value is 'atomic'.
    units : str 
        The LAMMPS units option associated with the data file. Default value
        is 'metal'.
    velocitiesstart : int or None
        The line of the file where the Velocities table content starts.

    Returns
    -------
    atomman.System
        The system with atom properties from the Velocities table assigned.
    """
    if velocitiesstart is not None:
        prop_info = velocities_prop_info(atom_style, units)
        system = amload('table', data, box=system.box, system=system,
                        prop_info=prop_info, skiprows=velocitiesstart,
                        nrows=system.natoms, comment='#', header=None)


    return system
