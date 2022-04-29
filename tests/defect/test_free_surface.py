import numpy as np
import pytest

import atomman as am
import atomman.unitconvert as uc
from atomman import Box, Atoms, System
from atomman.defect import FreeSurface, free_surface_basis
from atomman.tools.miller import vector_crystal_to_cartesian


@pytest.fixture
def bcc_Fe() -> System:
    a = uc.set_in_units(2.866, "angstrom")
    box = Box.cubic(a)
    atoms = Atoms(pos=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    system = System(atoms=atoms, box=box, scale=True, symbols="Fe")
    return system


@pytest.fixture
def anatase_TiO2() -> System:
    # mp-390
    a = uc.set_in_units(3.80271000, "angstrom")
    c = uc.set_in_units(9.74775200, "angstrom")
    box = Box.tetragonal(a, c)
    pos = [
        [0.00000000, 0.50000000, 0.25000000],
        [0.00000000, 0.00000000, 0.00000000],
        [0.50000000, 0.00000000, 0.75000000],
        [0.50000000, 0.50000000, 0.50000000],
        [0.00000000, 0.50000000, 0.45616300],
        [0.50000000, 0.50000000, 0.29383700],
        [0.00000000, 0.50000000, 0.04383700],
        [0.00000000, 0.00000000, 0.20616300],
        [0.50000000, 0.00000000, 0.95616300],
        [0.00000000, 0.00000000, 0.79383700],
        [0.50000000, 0.00000000, 0.54383700],
        [0.50000000, 0.50000000, 0.70616300],
    ]
    atype = [1] * 4 + [2] * 8
    atoms = Atoms(pos=pos, atype=atype)
    symbols = ['Ti'] * 4 + ['O'] * 8
    system = System(atoms=atoms, box=box, scale=True, symbols=symbols)
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


def test_free_surface_basis_rhombohedroal(alpha_Fe2O3):
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

    v1 = vector_crystal_to_cartesian(uvws[0], box)
    v2 = vector_crystal_to_cartesian(uvws[1], box)
    v3 = vector_crystal_to_cartesian(uvws[2], box)

    # check if v1 and v2 are on (hkl) plane
    np.testing.assert_almost_equal(np.dot(v1, planenormal), 0)
    np.testing.assert_almost_equal(np.dot(v2, planenormal), 0)

    # check right-handedness
    assert np.dot(np.cross(v1, v2), v3) > 0

    # check if v3 is out of (hkl) plane
    assert np.dot(v3, planenormal) > 0


def test_unique_shifts_binary(anatase_TiO2):
    # Example from Section 2.4 of W. Sun and G. Cedar, Surface Science, 617, 53-59 (2013)
    # The two possible shifts (terminations) are identical under a glide operation.
    pytest.importorskip("spglib")

    hkl = [1, 0, 0]
    surface = FreeSurface(hkl, ucell=anatase_TiO2)
    unique_shifts = surface.unique_shifts()
    assert len(surface.shifts) == 2
    assert len(unique_shifts) == 1


CUBIC_HKL_LISTS = [
    [1, 0, 0],
    [1, 1, 0],
    [1, 1, 1],
    [2, 1, 0],
    [2, 1, 1],
    [3, 1, 0],
    [3, 1, 1],
    [3, 2, 0],
    [3, 2, 1],
    [3, 2, 2],
    [3, 3, 1],
    [3, 3, 2],
]
TETRAGONAL_HKL_LISTS = [
    [0, 0, 1],
    [1, 0, 0],
    [1, 0, 1],
    [1, 1, 0],
    [1, 1, 1],
    [1, 1, 2],
    [2, 0, 1],
    [2, 1, 1],
]
HEXGONAL_HKL_LISTS = [
    [0, 0, 0, 1],
    [1, 0, -1, 0],
    [1, 0, -1, 1],
    [1, 0, -1, 2],
    [1, 1, -2, 0],
    [1, 1, -2, 1],
    [2, 0, -2, 1],
    [2, -1, -1, 2],
    [2, 1, -3, 0],
    [2, 1, -3, 1],
    [2, 1, -3, 2],
    [2, 2, -4, 1],
]

@pytest.mark.parametrize(
    "prototype,list_hkl,list_expected", [
        ('A1--Cu--fcc', CUBIC_HKL_LISTS, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        ('A2--W--bcc', CUBIC_HKL_LISTS, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        ("A3--Mg--hcp", HEXGONAL_HKL_LISTS, [1, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 1]),
        ("A4--C--dc", CUBIC_HKL_LISTS, [1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1]),
        ("A5--beta-Sn", TETRAGONAL_HKL_LISTS, [1, 1, 2, 1, 1, 1, 2, 2]),
        ('A6--In--bct', TETRAGONAL_HKL_LISTS, [1, 1, 1, 1, 1, 1, 1, 1]),
        ('A7--alpha-As', HEXGONAL_HKL_LISTS, [2, 1, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2]),
        ('Ah--alpha-Po--sc', CUBIC_HKL_LISTS, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        # TODO
        # ('A15--beta-W', CUBIC_HKL_LISTS, [2, 4, 2, 1, 2, 1, 1, 1, 2, 1, 1, 2]),
        # ("A3'--alpha-La--double-hcp", HEXGONAL_HKL_LISTS, [2, 3, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1]),
])
def test_unique_shifts_prototype(prototype, list_hkl, list_expected):
    # For some prototypes, the number of detected free surfaces differs from https://github.com/lmhale99/potentials-library/tree/master/free_surface.
    # A3'(alpha-La) (0001):
    #   seems to have two free surfaces, La(2a) surface and La(2c) one. but equivalent as slab models
    # - A4(diamond) (322) seems to have two free surfaces.
    # - A3' (100) seems to have two free surfaces.

    assert len(list_hkl) == len(list_expected)
    ucell = am.load("prototype", id=prototype)
    for hkl, expected in zip(list_hkl, list_expected):
        surface = FreeSurface(hkl, ucell)
        unique_shifts = surface.unique_shifts()
        assert len(unique_shifts) == expected
