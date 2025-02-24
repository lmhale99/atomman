import pytest
import numpy as np
import warnings

from atomman.plot.interpolate_contour import __grid_interpolate


def test_grid_interpolate():
    x = np.array([0, 0, 1])
    y = np.array([0, 1, 1])
    v = np.array([0, 1, 2])
    xlim = (0, 1)
    ylim = (0, 1)

    # extrapolation occurs
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid = __grid_interpolate(x, y, v, xlim=xlim, ylim=ylim)
    assert np.any(np.isnan(grid))

    # fill in extrapolated points
    fill_value = 1234
    filled_grid = __grid_interpolate(x, y, v, xlim=xlim, ylim=ylim, fill_value=fill_value)
    np.testing.assert_almost_equal(filled_grid[np.isnan(grid)], fill_value)
