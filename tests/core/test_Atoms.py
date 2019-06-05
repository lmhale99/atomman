# coding: utf-8

# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am

class Test_Atoms:
    def test_default(self):
        atoms = am.Atoms()
        assert atoms.natoms == 1
        assert atoms.natypes == 1
        assert atoms.atype[0] == 1
        assert np.allclose(atoms.pos[0], np.zeros(3))
        assert len(atoms) == atoms.natoms

    def build_example(self):
        atoms = am.Atoms(pos=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], atype=[2,1])
        atoms.test1 = np.array(['a', 'b'])
        atoms.prop(key='test2', value=[np.zeros((3,3))])
        atoms.view['test3'] = [1, 1]
        atoms.prop_atype('charge', [-1, 1])
        return atoms

    def test_assignments(self):
        atoms = self.build_example()
        assert atoms.natoms == 2
        assert atoms.natypes == 2
        assert atoms.atypes[0] == 1
        assert atoms.atypes[1] == 2
        assert atoms.prop(key='test1', index=0) == 'a'
        assert atoms.test1[1] == 'b'
        assert np.allclose(atoms.view['test2'][0], np.zeros((3,3)))
        assert np.allclose(atoms.test2[1], np.zeros((3,3)))
        assert atoms.view['test3'][0] == 1
        assert atoms.test3[1] == 1
        assert atoms.charge[0] == 1
        assert atoms.charge[1] == -1

    def test_df(self):
        atoms = self.build_example()
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
        assert 'test3' in df
        assert 'charge' in df

    def test_model(self):
        atoms = self.build_example()
        model = atoms.model(prop_unit={'atype':None, 'pos':'scaled',
                            'test1':None, 'test2':None, 'test3':'J',
                            'charge':'GPa'})
        atoms2 = am.Atoms(model=model)
        for prop in ['atype', 'pos', 'test2', 'test3', 'charge']:
            assert np.allclose(atoms.view[prop], atoms2.view[prop])
        assert atoms.test1[0] == atoms2.test1[0]
        assert atoms.test1[1] == atoms2.test1[1]
