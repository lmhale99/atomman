# coding: utf-8

# Standard Python libraries
from typing import Optional, Tuple, Union

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

# atomman imports
import atomman.unitconvert as uc
from .. import System

def get_prop_values(system: System,
                    prop_name: str,
                    prop_index: Union[int, tuple, None] = None,
                    prop_magnitude: bool = False,
                    prop: Optional[npt.ArrayLike] = None,
                    prop_unit: Optional[str] = None,
                    ) -> Tuple[str, np.ndarray]:
    """
    Utility method to manage the prop parameters used by the various plotting
    tools.
    
    Parameters
    ----------
    system : atomman.System
        The system with the per-atom property that you want to plot.
    prop_name : str
        The name of the per-atom property that you want to plot.
    prop_index : int or tuple, optional
        Specifies which component of a multidimensional property to plot.  Not
        needed if the property is scalar.
    prop_magnitude : bool, optional
        If True, plots the per-atom magnitude of a vector property.  Cannot be
        combined with index.  Default value is False.
    prop : array-like object, optional
        Values for the per-atom property to plot.  If not given, values will
        be taken as the "name" property of system.
    prop_unit : str or None, optional
        The units to use for the property value being plotted.  Default value
        is None, in which no unit conversion is applied.
        
    Returns
    -------
    prop_name : str
        The updated property name.
    prop : np.ndarray
        The per-atom property values.
    """

    # Get property from system
    if prop is None:
        prop = system.atoms.view[prop_name]

    # Handle index
    if prop_index is not None:
        if prop_magnitude is True:
            raise ValueError('prop_index and prop_magnitude cannot be combined')
        if isinstance(prop_index, (int, np.integer)):
            prop_index = [prop_index]
        else:
            prop_index = list(prop_index)
        for i in prop_index:
            prop_name += f'[{i+1}]'
        prop = prop[tuple([Ellipsis] + prop_index)]

    # Handle magnitude
    elif prop_magnitude is True:
        prop = np.linalg.norm(prop, axis=1)
        prop_name += '_mag'

    # Check that the property is of the right shape
    if prop.shape != (system.natoms, ):
        raise ValueError('property to plot must be a per-atom scalar')

    # Convert to the specified units
    prop = uc.get_in_units(prop, prop_unit)

    return prop_name, prop
