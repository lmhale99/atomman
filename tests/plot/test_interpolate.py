import pytest
import numpy as np

from atomman.plot.interpolate_contour import grid_interpolate_2d


def test_grid_interpolate_2d():
    x = np.array([0, 0, 1])
    y = np.array([0, 1, 1])
    v = np.array([0, 1, 2])
    range = [[0, 1], [0, 1]]

    # extrapolation occurs
    grid, _, _ = grid_interpolate_2d(x, y, v, range=range)
    assert np.any(np.isnan(grid))

    # fill in extrapolated points
    fill_value = 1234
    filled_grid, _, _ = grid_interpolate_2d(x, y, v, range=range, fill_value=fill_value)
    np.testing.assert_almost_equal(filled_grid[np.isnan(grid)], fill_value)
