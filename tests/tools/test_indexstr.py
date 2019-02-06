import pytest
import atomman as am 

def compile_indices(shape):
    alli = []
    allistr = []
    for i, istr in am.tools.indexstr(shape):
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
    