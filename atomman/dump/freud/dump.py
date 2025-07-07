from typing import Tuple

import freud

import numpy as np
import numpy.typing as npt

from ... import System

def dump(system) -> Tuple[freud.box.Box, np.ndarray]:
    """
    Convert an atomman.system into a freud.box.Box and an array of coordinate points. 

    NOTE: The generated freud representation should only be used for 
    interfacing with freud methods as it is an incomplete representation of the
    atomic system that ignores all per-atom properties except pos, a non-zero
    box origin, and the atomic positions are dropped to 32-bit precision.

    Returns
    -------
    box : freud.box.Box
        The freud representation of the cell box.
    points : numpy.ndarray
        The relative atomic positions expressed as an (N, 3) array of float32
        values.
    """
    # Convert atomman's box to freud's box
    lx = system.box.lx
    ly = system.box.ly
    lz = system.box.lz
    xy = system.box.xy / ly
    xz = system.box.xz / lz
    yz = system.box.yz / lz
    box = freud.box.Box(lx, ly, lz, xy, xz, yz)
    box.periodic = system.pbc.tolist()

    # Convert atomman's atom pos to freud's points
    points = system.atoms_prop('pos', scale=True)#.astype('float32')
    points = box.make_absolute(points)
    return (box, points)

