import datetime
from typing import Optional

from ..lammps_lib_parameters import dump as dump_lammps_lib_parameters
from ...typing import lammps

def dump(
    system,
    lmp: lammps.lammps,
    potential,
    prompt: bool = False,
    comments: bool = True,
    lammps_date: Optional[datetime.date] = None,
    tilt_large: bool = False,
    include_velocities: bool = True,
    no_atoms: bool = False,
    clear: bool = True):
    """
    Extracts data stored in an atomman System and passes it to a dynamic
    LAMMPS object to define basic settings, box, atoms, and potential
    information.
    
    Parameters
    ----------
    system : atomman.System 
        The system that LAMMPS will replicate.
    lmp :  lammps.lammps
        The LAMMPS object to which the system information will be passed.
    potential : atomman.lammps.Potential
        The potential object to use with the system.
    prompt : bool, optional
        If True, then a screen prompt will appear for radioactive elements
        with no standard mass to ask for the isotope to use. If False
        (default), then the most stable isotope will be automatically used.
    comments : bool, optional
        Indicates if print command lines detailing information on the potential
        are to be included.  Default value is True.
    lammps_date : datetime, optional
        The LAMMPS version date that is to be used.  If not given, will check
        the version date for the lmp object.
    tilt_large : bool, optional
        If True, then the box command will include the tilt large option for
        older LAMMPS versions to allow for larger tilt factors.  Default value
        is False.
    include_velocities : bool, optional
        Indicates if velocity information in the system (if present) is to
        be extracted and included in the returned outputs.  Default value is
        True.
    no_atoms : bool, optional
        Setting this to True will exclude atom definitions from the generated
        LAMMPS commands.  This allows for a box and potential to be set up
        with atoms being defined later.
    """
    # Extract the system parameters
    system_params = dump_lammps_lib_parameters(system, units=potential.units, 
                                               include_velocities=include_velocities)
    
    if no_atoms:
        system_params['natoms'] = 0

    # Pass the system parameters to the potential and run the associated LAMMPS commands
    potential.pair_lib_info(lmp, prompt=prompt, comments=comments,
                            lammps_date=lammps_date, tilt_large=tilt_large,
                            clear=clear, **system_params)
