import numpy as np

from atomman.region import Plane


def test_isclose():
    plane1 = Plane(normal=[1, 0, 0], point=[1, 0, 0])
    plane2 = Plane(normal=[1.0001, 0, 0], point=[1.0001, 1, 1])
    assert plane1 != plane2
    assert plane1.isclose(plane2, rtol=1e-4, atol=1e-4)
    assert plane2.isclose(plane1, rtol=1e-4, atol=1e-4)


def test_equality():
    plane1 = Plane(normal=[1, 0, 0], point=[1, 0, 0])
    plane2 = Plane(normal=[1, 0, 0], point=[1, 1, 1])
    plane3 = Plane(normal=[-1, 0, 0], point=[1, 0, 0])
    assert plane1 == plane2
    assert plane1 != plane3


def test_operate():
    plane = Plane(normal=[1, 0, 0], point=[1, 0, 0])
    # Rotate 90 degree along z-axis
    rotation = [
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1],
    ]
    translation = [0, 0, 0.5]
    new_plane = plane.operate(rotation, translation)
    np.testing.assert_almost_equal(new_plane.normal, [0, 1, 0])
    np.testing.assert_almost_equal(new_plane.point, [0, 1, 0.5])
