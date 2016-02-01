import numpy as np

def axes_check(axes, tol=1e-8):
    """Checks that axes are orthogonal and right handed, and returns the unit vectors."""
    uaxes = (axes.T / np.linalg.norm(axes, axis=1)).T
    
    assert np.allclose(np.dot(uaxes, uaxes.T), np.identity(3), atol=tol), 'axes are not orthogonal'
    assert np.allclose(np.cross(uaxes[0], uaxes[1]), uaxes[2], atol=tol), 'axes are not right handed'
    
    return uaxes