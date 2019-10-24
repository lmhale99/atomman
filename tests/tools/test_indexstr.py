# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

from atomman.tools import indexstr

def compile_indices(shape):
    alli = []
    allistr = []
    for i, istr in indexstr(shape):
        alli.append(i)
        allistr.append(istr)
    return alli, allistr

def test_indexstr():
    shape = (0,)
    alli, allistr = compile_indices(shape)
    assert alli == []
    assert allistr == []
    
    shape = (2,2,1)
    alli, allistr = compile_indices(shape)
    assert alli == [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)]
    assert allistr == ['[0][0][0]', '[0][1][0]', '[1][0][0]', '[1][1][0]']
    