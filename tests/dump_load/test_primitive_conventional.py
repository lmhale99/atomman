# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am

class Test_primitive_conventional:

    def conventional_ucell_builder(self, family, setting):
        """
        Construct a test unit cell based on family and setting

        Parameters
        ----------
        family :  str
            The crystal family: m, o, t, c
        setting : str
            The conventional basis setting: p, i, f, a, b, c
        """

        box = {}
        box['m'] = am.Box.monoclinic(a=3, b=4, c=5, beta=112)
        box['o'] = am.Box.orthorhombic(a=3, b=4, c=5)
        box['t'] = am.Box.tetragonal(a=3, c=5)
        box['c'] = am.Box.cubic(a=4)

        atoms = {}
        atoms['p'] = am.Atoms(pos=[[0.0, 0.0, 0.0]])
        atoms['a'] = am.Atoms(pos=[[0.0, 0.0, 0.0],
                                   [0.0, 0.5, 0.5]])
        atoms['b'] = am.Atoms(pos=[[0.0, 0.0, 0.0],
                                   [0.5, 0.0, 0.5]])
        atoms['c'] = am.Atoms(pos=[[0.0, 0.0, 0.0],
                                   [0.5, 0.5, 0.0]])
        atoms['i'] = am.Atoms(pos=[[0.0, 0.0, 0.0],
                                   [0.5, 0.5, 0.5]])
        atoms['f'] = am.Atoms(pos=[[0.0, 0.0, 0.0],
                                   [0.0, 0.5, 0.5],
                                   [0.5, 0.0, 0.5],
                                   [0.5, 0.5, 0.0]])

        return am.System(atoms=atoms[setting], box=box[family], symbols='Au',
                         scale=True, safecopy=True)

    def transform_and_validate(self, family, setting):

        # Create conventional cell and convert to and from primitive cell
        c_ucell1 = self.conventional_ucell_builder(family, setting)
        p_ucell = c_ucell1.dump('conventional_to_primitive', setting=setting)
        c_ucell2 = p_ucell.dump('primitive_to_conventional', setting=setting)

        # Check volume per atom
        assert np.allclose(c_ucell1.box.volume / c_ucell1.natoms,
                           p_ucell.box.volume / p_ucell.natoms) 
        
        # Check box and atom positions of the two conventional cells
        assert np.allclose(c_ucell1.box.vects, c_ucell2.box.vects)
        aindex1, aindex2 = np.meshgrid(range(c_ucell1.natoms), range(c_ucell1.natoms))
        assert np.isclose(c_ucell1.dmag(aindex1.flatten(), aindex2.flatten()), 0.0).sum() == c_ucell1.natoms

    def test_bcc(self):
        self.transform_and_validate('c', 'i')

    def test_fcc(self):
        self.transform_and_validate('c', 'f')

    def test_bct(self):
        self.transform_and_validate('t', 'i')

    def test_fco(self):
        self.transform_and_validate('o', 'f')

    def test_bco(self):
        self.transform_and_validate('o', 'i')

    def test_sco(self):
        self.transform_and_validate('o', 'a')
        self.transform_and_validate('o', 'b')
        self.transform_and_validate('o', 'c')

    def test_scm(self):
        self.transform_and_validate('m', 'a')
        self.transform_and_validate('m', 'b')
        self.transform_and_validate('m', 'c')