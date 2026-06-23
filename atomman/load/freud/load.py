# Standard Python imports
from typing import Optional

# http://www.numpy.org/
import numpy.typing as npt

# https://freud.readthedocs.io/en/latest/index.html
import freud

# atomman imports
from ... import Atoms, Box, System

def load(freud_box: freud.box.Box,
         positions: npt.ArrayLike,
         origin: Optional[npt.ArrayLike] = None,
         atype: Optional[npt.ArrayLike] = None,
         symbols: Optional[tuple] = None) -> System:
    """
    Converts the freud box, positions representation into an atomman System object.

    Parameters
    ----------
    box : freud.box.Box
        A cell box defined as a freud Box object.
    positions : array-like
        (Nx3) array of positions given in the freud format.
    origin : array-like or None, optional
        An origin for the atomman Box representation to position it in absolute
        space.  If not given, the origin will be taken as [0,0,0].
    atype : array-like or None, optional
        Array of length N of integer atom types to associate with each position.
        If not given, all atype values will be 1. 
    symbols: tuple or None, optional
        Optional elemental symbols to assign for each unique atype value.
    """

    # Convert the freud Box to an atomman Box
    lx = freud_box.Lx
    ly = freud_box.Ly
    lz = freud_box.Lz
    xy = freud_box.xy * ly
    xz = freud_box.xz * lz
    yz = freud_box.yz * lz
    box = Box(lx=lx, ly=ly, lz=lz, xy=xy, xz=xz, yz=yz, origin=origin)
    
    # Build atomman Atoms from freud positions (and atype)
    pos = freud_box.make_fractional(positions)
    atoms = Atoms(atype=atype, pos=pos)

    # Build system and assign pbc and symbols
    system = System(atoms=atoms, box=box, symbols=symbols,
                    pbc=freud_box.periodic, scale=True)

    return system
    