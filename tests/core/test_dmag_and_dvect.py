# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am

def test_dmag_and_dvect():
    pos_0 = np.random.uniform(1.0, 10.0, size=(200,3))
    pos_0[10] = [12.1234, 12.1234, 12.1234]
    pos_0[11] = [7.124, 132.15, 12.235]

    pos_1 = np.random.uniform(1.0, 10.0, size=(200,3))
    pos_1[10] = [12.1234, 12.1234, 12.1234]
    pos_1[11] = [18.124, 1.15, 1.235]

    box = am.Box(vects=[[ 20.0,   0.0,   0.0],
                        [  0.0, 135.0,   0.0],
                        [  0.0,   0.0,  20.0]])

    # Test for all periodic boundaries
    pbc = [True, True, True]
    dvect = am.dvect(pos_0, pos_1, box, pbc)
    dmag = am.dmag(pos_0, pos_1, box, pbc)

    assert np.allclose(np.linalg.norm(dvect, axis=1), dmag)
    assert np.allclose(dvect[10], [0.0, 0.0, 0.0])
    assert np.allclose(dvect[11], [-9.,  4.,  9.])
    assert np.isclose(dmag[10], 0.0)
    assert np.isclose(dmag[11], 13.341664064126334)

    # Test for all nono-periodic boundaries
    pbc = [False, False, False]
    dvect = am.dvect(pos_0, pos_1, box, pbc)
    dmag = am.dmag(pos_0, pos_1, box, pbc)

    assert np.allclose(np.linalg.norm(dvect, axis=1), dmag)
    assert np.allclose(dvect[10], [0.0, 0.0, 0.0])
    assert np.allclose(dvect[11], [11., -131.,  -11.])
    assert np.isclose(dmag[10], 0.0)
    assert np.isclose(dmag[11], 131.92043056327552)