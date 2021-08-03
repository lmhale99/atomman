# coding: utf-8

from DataModelDict import DataModelDict as DM

from potentials.record import PotentialLAMMPS, PotentialLAMMPSKIM

def Potential(model, name=None, pot_dir=None, kim_id=None):
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
    
    Returns
    -------
    potentials.record.PotentialLAMMPS
        Returned if model contains potential-LAMMPS content.
    potentials.record.PotentialLAMMPSKIM
        Returned if model contains potential-LAMMPS-KIM content.
    """
    model = DM(model)
    try:
        # Search for potential-LAMMPS branch
        model.find('potential-LAMMPS')
    except:
        try:
            # Search for potential-LAMMPS-KIM branch
            model.find('potential-LAMMPS-KIM')
        except:
            raise ValueError('Failed to find either potential-LAMMPS or potential-LAMMPS-KIM content')
        else:
            return PotentialLAMMPSKIM(model=model, name=name, id=kim_id)
    else:
        return PotentialLAMMPS(model=model, name=name, pot_dir=pot_dir)