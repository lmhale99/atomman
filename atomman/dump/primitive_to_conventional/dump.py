# coding: utf-8
# Standard Python libraries
from typing import Union, Tuple

# http://www.numpy.org/
import numpy as np

# atomman imports
from ... import System
from ...tools import miller

def dump(system: System,
         setting: str = 'p',
         return_transform: bool = False
         ) -> Union[System, Tuple[System, np.ndarray]]:
    """
    Transforms a primitive unit cell system into a conventional unit cell
    system of the given Bravais space lattice setting.  The
    primitive_to_conventional and conventional_to_primitive dump styles are
    meant to be inverse operations, within floating point precision, to provide
    compatible primitive and conventional unit cells.

    NOTE: This dump style expects that the original starting system is a
    primitive unit cell, although no checks are performed to assert this!

    Parameters
    ----------
    system : atomman.System
        A primitive unit cell system to find the corresponding conventional
        unit cell for.
    setting : str, optional
        The conventional cell space group lattice setting. Allowed values are
        'p' for primitive, 'f' for face-centered, 'i' for body centered, and
        'a', 'b', or 'c' for side-centered.
    return_transform : bool, optional
        Indicates if the Cartesian transformation matrix associated with
        rotating from the primitive cell to conventional cell orientations is
        returned.  Default value is False.

    Returns
    -------
    c_ucell : atomman.System
        The conventional unit cell obtained by transforming the given primitive
        unit cell system.
    transform : numpy.ndarray
        The Cartesian transformation matrix associated with converting from the
        primitive cell orientation to the conventional cell orientation.  Only
        returned if return_transform is True.
    """

    # Get rotations to convert primitive to conventional
    p2c_uvws = miller.vector_conventional_to_primitive(np.identity(3),
                                                       setting=setting)

    c_ucell, transform = system.rotate(p2c_uvws, return_transform=True)
    c_ucell.wrap()

    if return_transform:
        return c_ucell, transform
    else:
        return c_ucell
