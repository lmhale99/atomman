from typing import Optional
import secrets
              
seedmax: int = 900000000
"""int: The max random number generator seed LAMMPS supports"""

def newseed() -> int:
    """Randomly selects a random number generator seed for LAMMPS"""
    # A LAMMPS seed must be > 0
    return secrets.randbelow(seedmax - 1) + 1

def seed(val: Optional[int] = None) -> int:
    """
    Either generates a new random number generator seed for LAMMPS or checks if
    a given int is a valid seed number for LAMMPS.

    Parameters
    ----------
    val : int or None, optional
        If int, checks that it is in the proper range.  If None (default), a
        new seed will be randomly selected using secrets.
    
    Returns
    -------
    seed : int
        Same as val if it was an acceptable int, or a new seed if val was None.
    """
    if val is None:
        # Generate new seed
        return newseed()

    # Convert to int
    val = int(val)
    
    # Check that val is in the acceptable range
    if val <= 0 :
        raise ValueError('seed value must be > 0')
    elif val > seedmax:
        raise ValueError(f'seed value must be <= {seedmax}')
    
    return val