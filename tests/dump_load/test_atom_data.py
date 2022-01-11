# coding: utf-8

# https://docs.pytest.org/en/latest/
import pytest

# http://www.numpy.org/
import numpy as np

import atomman as am
from atomman.load import FileFormatError

class Test_atom_data:

    def load_dump(self, system1, atom_style=None, units=None, potential=None):
        """
        Utility function that dumps, loads, and dumps, then asserts two dumps
        are equivalent
        """
        # dump system1 to content1
        content1 = system1.dump('atom_data', atom_style=atom_style,
                                units=units, potential=potential,
                                return_info=False)

        # load content1 to system2
        system2 = am.load('atom_data', content1, pbc=system1.pbc,
                          symbols=system1.symbols, atom_style=atom_style,
                          units=units)

        # dump system2 to content2
        content2 = system2.dump('atom_data', atom_style=atom_style,
                                units=units, potential=potential,
                                return_info=False)

        # Check that the two dumps are equivalent
        assert content1 == content2

    def test_badfilename(self):
        """Raise FileNotFoundError for missing file"""
        with pytest.raises(FileNotFoundError):
            am.load('atom_data', 'nofile.dat')

    @property
    def data_lines(self):
        """List of file lines for testing"""
        return ['',
                '4 atoms',
                '2 atom types',
                '0.0000000000000 1.2569401286715 xlo xhi',
                '0.0000000000000 4.0169045159655 ylo yhi',
                '0.0000000000000 0.7639194529387 zlo zhi',
                '1.1755124391923 0.2312234378152 0.0768582022466 xy xz yz',
                '',
                'Atoms',
                '',
                '1 1 1.5543211985914 2.8647116844039 0.1665374455367',
                '2 1 2.3982049184958 4.0169567962790 0.4422255026206',
                '3 1 1.6332403214635 1.8356166738031 0.3819388814516',
                '4 1 1.2352186803268 0.3081624624139 0.5885979440072']

    def test_goodfile(self):
        """Test that full data_lines is valid"""
        content = '\n'.join(self.data_lines)
        am.load('atom_data', content)

    def test_badfile_no_natoms(self):
        """Raise FileFormatError if # atoms line is missing"""

        content = '\n'.join(self.data_lines[:1] + self.data_lines[2:])
        with pytest.raises(FileFormatError):
            am.load('atom_data', content)

    def test_badfile_no_box(self):
        """Raise FileFormatError if any *lo *hl lines are missing"""

        # missing xlo xhi
        content = '\n'.join(self.data_lines[:3] + self.data_lines[4:])
        with pytest.raises(FileFormatError):
            am.load('atom_data', content)

        # missing ylo yhi
        content = '\n'.join(self.data_lines[:4] + self.data_lines[5:])
        with pytest.raises(FileFormatError):
            am.load('atom_data', content)

        # missing zlo zhi
        content = '\n'.join(self.data_lines[:5] + self.data_lines[6:])
        with pytest.raises(FileFormatError):
            am.load('atom_data', content)

    def test_badfile_no_atoms(self):
        """Raise FileFormatError if Atoms table is missing"""

        content = '\n'.join(self.data_lines[:8])
        with pytest.raises(FileFormatError):
            am.load('atom_data', content)

    def test_atomic_no_imageflags(self):

        box = am.Box(vects=[[1.25694013, 0.        , 0.        ],
                            [1.17551244, 4.01690452, 0.        ],
                            [0.23122344, 0.0768582 , 0.76391945]])
        pos = [[0.53342537, 0.70899278, 0.21800393],
               [0.8766086 , 0.98893671, 0.57889022],
               [0.78898178, 0.44740662, 0.49997271],
               [0.78302097, 0.06197394, 0.77049739]]
        atoms = am.Atoms(pos=pos, atype=1)
        symbols = 'Al'
        pbc = (True, True, True)

        system = am.System(atoms=atoms, box=box, scale=True,
                           symbols=symbols, pbc=pbc)

        self.load_dump(system, atom_style='atomic', units='metal')

    def test_atomic_imageflags(self):

        box = am.Box(vects=[[1.25694013, 0.        , 0.        ],
                            [1.17551244, 4.01690452, 0.        ],
                            [0.23122344, 0.0768582 , 0.76391945]])
        pos = [[1.53342537, 0.70899278, 0.21800393],
               [0.8766086 , 1.98893671, 0.57889022],
               [0.78898178, 0.44740662, 1.49997271],
               [-2.78302097, 1.06197394, 0.77049739]]
        atoms = am.Atoms(pos=pos, atype=1)
        symbols = 'Al'
        pbc = (True, True, True)

        system = am.System(atoms=atoms, box=box, scale=True,
                           symbols=symbols, pbc=pbc)

        self.load_dump(system, atom_style='atomic', units='metal')
