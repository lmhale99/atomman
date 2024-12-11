# coding: utf-8

# Standard Python libraries
import io
from typing import Optional, Union

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/potentials
from potentials import load_record
from potentials.record.PotentialLAMMPS import PotentialLAMMPS
from potentials.record.PotentialLAMMPSKIM import PotentialLAMMPSKIM


def Potential(model: Union[str, io.IOBase],
              name: Optional[str] = None,
              pot_dir: Optional[str] = None,
              kim_id: Optional[str] = None,
              potkey: Optional[str] = None,
              potid: Optional[str] = None,
              symbolset: Union[str, list, None] = None
              ):
    """
    Class initializer switch for PotentialLAMMPS or PotentialLAMMPSKIM objects.

    Parameters
    ----------
    model : str or file-like object
        A JSON/XML data model containing a potential-LAMMPS or a
        potential-LAMMPS-KIM branch.
    name : str, optional
        The record name to use.  If not given, this will be set to the
        potential's id.
    pot_dir : str, optional
        The path to a directory containing any artifacts associated with
        the potential.  Default value is None, which assumes any required
        files will be in the working directory when LAMMPS is executed.
    kim_id : str, optional
        The full KIM model id indicating the version of a KIM model to use.
        If not given, then the newest known version will be used.
    potkey : str, optional
        Specifies which potential (by potkey value) to use.  Only important
        if it is a kim model associated with multiple potential entries.
        This controls which symbols are available.
    potid : str, optional
        Specifies which potential (by potid value) to use.  Only important
        if it is a kim model associated with multiple potential entries.
        This controls which symbols are available.
    symbolset : str or list, optional
        Specifies which potential (by symbols value) to use.  Only important
        if it is a kim model associated with multiple potential entries.
        If potkey or potid is not given, then the first potential entry
        found with all listed symbols will be selected.
    
    Returns
    -------
    potentials.record.BasePotentialLAMMPS.BasePotentialLAMMPS
        The matching LAMMPS Potential object
    """
    model = DM(model)
    try:
        # Search for potential-LAMMPS branch
        model.find('potential-LAMMPS')
    except ValueError:
        try:
            # Search for potential-LAMMPS-KIM branch
            model.find('potential-LAMMPS-KIM')
        except ValueError as err:
            raise ValueError('Failed to find either potential-LAMMPS or potential-LAMMPS-KIM content') from err
        else:
            return load_record('potential_LAMMPS_KIM', model=model, name=name, id=kim_id,
                               potkey=potkey, potid=potid, symbolset=symbolset)
    else:
        return load_record('potential_LAMMPS', model=model, name=name, pot_dir=pot_dir)