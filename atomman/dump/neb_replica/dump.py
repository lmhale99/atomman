# coding: utf-8

# Standard Python libraries
import io
from typing import Optional, Union

def dump(system,
         f: Union[str, io.IOBase, None] = None,
         id_key: Optional[str] = None,
         id_start0: bool = True,
         float_format: str = '%.13f') -> Optional[str]:
    """
    Generates the atomic coordinate file used by LAMMPS NEB for specifying
    intermediate and final replica configurations.

    NOTE that id_key SHOULD BE SET if the replica only includes the moved
    atoms subset of the reference initial/first replica system!!!

    Parameters
    ----------
    system : atomman.System
        The system to write to the NEB replica file.  
    f : str or file-like object, optional
        File path or file-like object to write the content to.  If not given,
        then the content is returned as a str.
    id_key : str, optional
        The name of the atoms property of system to use for the atomic ids.
        If not given, then the indices of the atoms will be used which should
        only be done if the atoms in the given system have the same number and
        order as the reference system used for the initial/first replica.  This
        atoms property should map the moved atoms to the corresponding ids of
        the initial/first replica. 
    id_start0 : bool, optional
        LAMMPS ids start at 1 whereas atomman uses atom indices which start at
        0.  If idstart0 is True (default) then this indicates that the id_key
        values are relative to the atomman atoms indices and should be
        increased by 1 when dumped.  If False, then the id_key values are used
        as is and assumed to be relative to the LAMMPS atom ids.
    float_format : str, optional
        c-style formatting string for floating point values.  Default value is
        '%.13f'.
    """
    # Build a dataframe for the system's atom info
    df = system.atoms_df()

    # Manage id_key values
    if id_key is None:
        # Create new 'atom_id' field
        df['atom_id'] = range(1, system.natoms+1)
        include_keys = ['atom_id', 'pos[0]', 'pos[1]', 'pos[2]']
    else:
        include_keys = [id_key, 'pos[0]', 'pos[1]', 'pos[2]']

        # Bump id field by 1 if needed
        if id_start0 is True:
            df[id_key] += 1
    
    # Generate the content
    content = f'{system.natoms}\n'
    content += df[include_keys].to_csv(sep=' ', index=None, header=False,
                                       float_format=float_format,
                                       lineterminator='\n')

    # Save to the file-like object
    if hasattr(f, 'write'):
        f.write(content)
    
    # Save to the file name
    elif f is not None:
        with open(f, 'w', encoding='UTF-8') as fp:
            fp.write(content)
    
    # Return as a string
    else:
        return content