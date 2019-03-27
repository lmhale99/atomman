import pytest
import atomman as am

# http://www.numpy.org/
import numpy as np

class Test_Atoms:
    def test_default(self):
        atoms = am.Atoms()
        assert atoms.natoms == 1
        assert atoms.natypes == 1
        assert atoms.atype[0] == 1
        assert np.allclose(atoms.pos[0], np.zeros(3))
        assert len(atoms) == atoms.natoms

    def test_1(self):
        atoms = am.Atoms(pos=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], atype=[2,1])
        assert atoms.natoms == 2
        assert atoms.natypes == 2
        assert atoms.atypes[0] == 1
        assert atoms.atypes[1] == 2
        
        atoms.test1 = np.array(['a', 'b'])
        assert atoms.prop(key='test1', index=0) == 'a'
        assert atoms.test1[1] == 'b'
        
        atoms.prop(key='test2', value=[np.zeros((3,3))])
        assert np.allclose(atoms.view['test2'][0], np.zeros((3,3)))
        assert np.allclose(atoms.test2[1], np.zeros((3,3)))

        atoms.view['test3'] = [1, 1]
        assert atoms.view['test3'][0] == 1
        assert atoms.test3[1] == 1
        
        atoms.prop_atype('charge', [-1, 1])       
        assert atoms.charge[0] == 1
        assert atoms.charge[1] == -1

        df = atoms.df()
        assert len(df) == 2
        assert 'atype' in df
        assert 'pos[0]' in df
        assert 'pos[1]' in df
        assert 'pos[2]' in df
        assert 'test1' in df
        assert 'test2[0][0]' in df
        assert 'test2[0][1]' in df
        assert 'test2[0][2]' in df
        assert 'test2[1][0]' in df
        assert 'test2[1][1]' in df
        assert 'test2[1][2]' in df
        assert 'test2[2][0]' in df
        assert 'test2[2][1]' in df
        assert 'test2[2][2]' in df
        assert 'charge' in df
