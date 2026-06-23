# Standard Python imports
from typing import Optional

# http://www.numpy.org/
import numpy as np

def random_rotation(x: Optional[float] = None,
                    y: Optional[float] = None,
                    z: Optional[float] = None,
                    rng: Optional[np.random.Generator] = None,
                    seed: Optional[int] = None):
    """
    Selects a random rotation matrix using a method that uniformly samples
    the phase space.  See Graphics Gems III (1992) p117-120 and p463-464.

    This implementation uses numpy to make the math in the code clearer,
    with a minor performance hit compared to the original "c optimized" code.

    To create small uniform "perturbations", you can limit x and z to a max value < 1.

    Parameters
    ----------
    x : float or None, optional
        Random number in range [0, 1] used to select the initial rotation term.
        If None (default), will select a value using numpy.random.rand().
    y : float or None, optional
        Random number in range [0, 1] used to select the pole deflection direction term.
        If None (default), will select a value using numpy.random.rand().
    z : float or None, optional
        Random number in range [0, 1] used to select the pole deflection magnitude term.
        If None (default), will select a value using numpy.random.rand().
    rng : numpy.random.Generator or None, optional
        An existing random number generator to use for the x,y,z selection.  If None
        (default), a new generator will be created based on seed if needed.
    seed : int or None, optional
        A random number seed to use when creating a generator if rng is not given.
        If None (default), then then fresh, unpredictable entropy will be pulled from the OS.
    """
    # Create a new random number generator if needed
    if rng is None:
        rng = np.random.default_rng(seed)

    # Select three random variables between 0 and 1 as needed
    if x is None:
        x = rng.random()
    if y is None:
        y = rng.random()
    if z is None:
        z = rng.random()

    # Scale the random variables
    theta = 2 * np.pi * x   # Rotation about original pole z
    phi = 2 * np.pi * y     # Direction of pole deflection
    z = 2 * z               # Magnitude of pole deflection
    
    # R is rotation matrix around pole z
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    R = np.array([
        [ cos_theta, sin_theta, 0.0],
        [-sin_theta, cos_theta, 0.0],
        [       0.0,       0.0, 1.0]])
    
    # V is vector for distributing points over the sphere
    # Note that it is scaled to cancel Householder Matrix stuff...
    V = np.array([np.sin(phi) * z**0.5, np.cos(phi) * z**0.5, (2.0 - z) ** 0.5])
   
    # S is row vector S = Transpose(V) * R
    S = np.dot(V, R)
    
    # Transformation matrix T = V S - R
    T = np.outer(V, S) - R

    return T