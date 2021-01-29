# coding: utf-8

# http://www.numpy.org/
import numpy as np

def central_difference(fxn, coord, shift=1e-5):
    """
    Computes the gradient of a function at a set of coordinates.
    
    Parameters
    ----------
    fxn : function
        The function to compute gradients for.
    coord : array-like object
        The coordinates to evaluate the gradient at. The number of
        derivatives calculated will be based on the size of the
        coordinate's final dimension, i.e. if coord is (10, 3), then three
        derivatives will be computed for each of the ten coordinates.
    shift : float
        The shift step size to use when evaluating the derivatives.
    
    Returns
    -------
    gradient : array-like object
        The gradient array with the same shape as coord.
    """
    coord = np.asarray(coord)
    
    # Identify number of derivative dimensions based on final array dimension
    ndim = coord.shape[-1]
    
    # Build gradient array
    gradient = np.zeros_like(coord)
    fps = []
    fms = []
    for i in range(ndim):
        
        # Create δ vector to displace along one dimension
        δ = np.zeros(ndim)
        δ[i] = shift
        
        # Compute gradients for that dimension
        fplus = fxn(coord + δ)
        fminus = fxn(coord - δ)
        fps.append(fplus)
        fms.append(fminus)
        
        gradient[..., i] = (fplus - fminus) / (2 * shift)
    
    return gradient