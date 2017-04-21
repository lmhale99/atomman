import numpy as np
from .supersize import supersize

def hex_to_ortho(cell, tol=1e-5):
    try:
        assert np.isclose(cell.box.a, cell.box.b, rtol=tol)
        assert not np.isclose(cell.box.a, cell.box.c, rtol=tol)
        assert np.isclose(cell.box.alpha, 90.0, rtol=tol)
        assert np.isclose(cell.box.beta, 90.0, rtol=tol)
        assert np.isclose(cell.box.gamma, 120.0, rtol=tol) or np.isclose(cell.box.gamma, 60.0, rtol=tol)
    except:
        raise ValueError('Not a standard hexagonal box')
    
    newcell = supersize(cell, 1, 2, 1) 
    newcell.normalize(style='lammps')
    
    new_b = cell.box.b * 3.**0.5
    newcell.box_set(origin=[0.0, 0.0, 0.0], scale=True)
    newcell.box_set(a=cell.box.a, b=new_b, c=cell.box.c, 
                    alpha=90.0, beta=90.0, gamma=90.0, 
                    origin=cell.box.origin)
    newcell.wrap()     
    return newcell