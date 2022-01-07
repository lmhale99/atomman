import numpy as np
import pytest

import atomman.unitconvert as uc
from atomman import Box, Atoms, System
from atomman.defect import FreeSurface, free_surface_basis


@pytest.fixture
def bcc_Fe() -> System:
    a = uc.set_in_units(2.866, "angstrom")
    box = Box.cubic(a)
    atoms = Atoms(pos=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    system = System(atoms=atoms, box=box, scale=True, symbols="Fe")
    return system


@pytest.fixture
def alpha_Fe2O3() -> System:
    """
    Working example from W. Sun and G. Cedar, Surface Science, 617, 53-59 (2013).
    """
    # mp-19770
    a = uc.get_in_units(5.10485225, "angstrom")
    c = uc.get_in_units(13.91329700, "angstrom")
    box = Box.hexagonal(a, c)

    pos = [
        [0.00000000, 0.00000000, 0.85417500],
        [0.00000000, 0.00000000, 0.14582500],
        [0.66666667, 0.33333333, 0.97915833],
        [0.33333333, 0.66666667, 0.02084167],
        [0.66666667, 0.33333333, 0.18750833],
        [0.66666667, 0.33333333, 0.47915833],
        [0.33333333, 0.66666667, 0.31249167],
        [0.00000000, 0.00000000, 0.35417500],
        [0.33333333, 0.66666667, 0.52084167],
        [0.33333333, 0.66666667, 0.81249167],
        [0.00000000, 0.00000000, 0.64582500],
        [0.66666667, 0.33333333, 0.68750833],
        [0.97154417, 0.33333333, 0.08333333],
        [0.02845583, 0.66666667, 0.91666667],
        [0.66666667, 0.63821083, 0.08333333],
        [0.33333333, 0.36178917, 0.91666667],
        [0.36178917, 0.02845583, 0.08333333],
        [0.63821083, 0.97154417, 0.91666667],
        [0.63821083, 0.66666667, 0.41666667],
        [0.69512250, 0.00000000, 0.25000000],
        [0.33333333, 0.97154417, 0.41666667],
        [0.00000000, 0.69512250, 0.25000000],
        [0.02845583, 0.36178917, 0.41666667],
        [0.30487750, 0.30487750, 0.25000000],
        [0.30487750, 0.00000000, 0.75000000],
        [0.36178917, 0.33333333, 0.58333333],
        [0.00000000, 0.30487750, 0.75000000],
        [0.66666667, 0.02845583, 0.58333333],
        [0.69512250, 0.69512250, 0.75000000],
        [0.97154417, 0.63821083, 0.58333333],
    ]
    atype = [1] * 12 + [2] * 18
    atoms = Atoms(atype=atype, pos=pos)

    symbols = ("Fe", "O")
    system = System(atoms=atoms, box=box, scale=True, symbols=symbols)
    return system


def test_free_surface_basis_cubic(bcc_Fe):
    box = bcc_Fe.box
    _test_free_surface_basis(box, [1, 0, 0])
    _test_free_surface_basis(box, [0, 2, 0])
    _test_free_surface_basis(box, [0, 0, 3])
    _test_free_surface_basis(box, [0, 1, -1])
    _test_free_surface_basis(box, [2, 0, 3])
    _test_free_surface_basis(box, [4, 5, 0])
    _test_free_surface_basis(box, [1, 2, -3])


def test_free_surface_basis_hexagonal(alpha_Fe2O3):
    box = alpha_Fe2O3.box
    _test_free_surface_basis(box, [1, 0, -1, 4])
    _test_free_surface_basis(box, [1, 0, 0])
    _test_free_surface_basis(box, [0, 2, 0])
    _test_free_surface_basis(box, [0, 0, 3])
    _test_free_surface_basis(box, [0, 1, -1])
    _test_free_surface_basis(box, [2, 0, 3])
    _test_free_surface_basis(box, [4, 5, 0])
    _test_free_surface_basis(box, [1, 2, -3])


def _test_free_surface_basis(box: Box, hkl):
    uvws, planenormal = free_surface_basis(
        hkl=hkl,
        box=box,
        return_hexagonal=False,
        return_planenormal=True,
    )

    v1 = box.cart(uvws[0])
    v2 = box.cart(uvws[1])
    v3 = box.cart(uvws[2])

    # check if v1 and v2 are on (hkl) plane
    np.testing.assert_almost_equal(np.dot(v1, planenormal), 0)
    np.testing.assert_almost_equal(np.dot(v2, planenormal), 0)

    # check right-handedness
    assert np.dot(np.cross(v1, v2), v3) > 0

    # check if v3 is out of (hkl) plane
    assert np.dot(v3, planenormal) > 0
