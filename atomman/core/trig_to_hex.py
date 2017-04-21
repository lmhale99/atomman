import numpy as np
from .supersize import supersize

def trig_to_hex(cell, tol=1e-5):    
    try:
        assert np.isclose(cell.box.a, cell.box.b, rtol=tol)
        assert np.isclose(cell.box.a, cell.box.c, rtol=tol)
        assert np.isclose(cell.box.alpha, cell.box.beta, rtol=tol)
        assert np.isclose(cell.box.alpha, cell.box.beta, rtol=tol)
        assert cell.box.alpha < 120.0
    except:
        raise ValueError('Not a standard trigonal box')
    
    newcell = supersize(cell, 3, 1, 1) 

    newcell.box_set(avect = cell.box.bvect - cell.box.avect,
                    bvect = cell.box.avect - cell.box.cvect,
                    cvect = cell.box.avect + cell.box.bvect + cell.box.cvect,
                    origin= cell.box.origin)
    newcell.wrap()   
    return newcell