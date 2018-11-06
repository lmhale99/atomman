import pytest
import atomman as am
import numpy as np

class Test_axes_check:
    def test_parallel(self):
        assert np.allclose(am.tools.axes_check(5*np.eye(3)), np.eye(3))
        
    def test_110(self):
        axes = np.array([[ 1, 1, 0], 
                         [-1, 1, 0], 
                         [ 0, 0, 1]])
        
        uaxes = np.array([[ 2.**0.5/2., 2.**0.5/2., 0], 
                          [-2.**0.5/2., 2.**0.5/2., 0], 
                          [ 0, 0, 1.0]])
        assert np.allclose(am.tools.axes_check(axes), uaxes)
        
    def test_bad_shape(self):
        with pytest.raises(AssertionError):
            am.tools.axes_check([[12,124], [124,62], [1234, 124]])
            
    def test_left_handed(self):
        axes = np.array([[ 1, 1, 0], 
                         [ 1,-1, 0], 
                         [ 0, 0, 1]])
        with pytest.raises(ValueError):
            am.tools.axes_check(axes)
            
    def test_not_orthogonal(self):
        axes = np.array([[ 1, 1, 0], 
                         [ 1,-3, 0], 
                         [ 0, 0, 1]])
        with pytest.raises(ValueError):
            am.tools.axes_check(axes)        