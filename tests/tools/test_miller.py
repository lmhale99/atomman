import pytest
import atomman as am 
import numpy as np 

def test_vector3to4():
    assert np.all(am.tools.miller.vector3to4([ 1, 2, 3]) == np.array([ 0, 1,-1, 3]))
    assert np.all(am.tools.miller.vector3to4([ 0, 1, 0]) == np.array([-1, 2,-1, 0]))
    assert np.all(am.tools.miller.vector3to4([-3, 0, 1]) == np.array([-2, 1, 1, 1]))
    assert np.all(am.tools.miller.vector3to4([-2, 2, 2]) == np.array([-1, 1, 0, 1]))
    with pytest.raises(AssertionError):
        am.tools.miller.vector3to4([ 1, 2, 3, 4])


def test_vector4to3():
    assert np.all(am.tools.miller.vector4to3([ 0, 1,-1, 3]) == np.array([ 1, 2, 3]))
    assert np.all(am.tools.miller.vector4to3([-1, 2,-1, 0]) == np.array([ 0, 1, 0]))
    assert np.all(am.tools.miller.vector4to3([-2, 1, 1, 1]) == np.array([-3, 0, 1]))
    assert np.all(am.tools.miller.vector4to3([-2, 2, 0, 2]) == np.array([-1, 1, 1]))
    with pytest.raises(AssertionError):
        am.tools.miller.vector4to3([ 1, 2, 3])

def test_vectortocartesian():
    # Test 3 indices
    box = am.Box(a=3, b=4, c=5, alpha=90, beta=90, gamma=90)
    assert np.allclose(am.tools.miller.vectortocartesian([1,2,3], box),
                       np.array([ 3, 8, 15]))

    box = am.Box(a=3, b=4, c=5, alpha=90, beta=90, gamma=80)
    assert np.allclose(am.tools.miller.vectortocartesian([1,1,1], box),
                       np.array([[3.69459271, 3.93923101, 5.        ]]))
    # Test 4 indices
    box = am.Box(a=3, b=3, c=5, alpha=90, beta=90, gamma=120)
    assert np.allclose(am.tools.miller.vectortocartesian([-2, 1, 1, 1], box),
                       np.array([-9, 0, 5]))
    
    box = am.Box(a=3, b=4, c=5, alpha=90, beta=90, gamma=120)
    with pytest.raises(AssertionError):
        am.tools.miller.vectortocartesian([-2, 1, 1, 1], box)