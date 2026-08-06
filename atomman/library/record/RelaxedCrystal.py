# coding: utf-8

# Standard Python imports
from typing import Tuple

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class RelaxedCrystal(Record):
    """
    Class for representing relaxed_crystal records that provide the structure
    information for crystal structures relaxed using a specific interatomic
    potential.
    """

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'relaxed_crystal'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'relaxed-crystal'

    @property
    def xsd_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsd schema"""
        return ('atomman.library.xsd', f'{self.style}.xsd')

    @property
    def xsl_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsl transformer"""
        return ('atomman.library.xsl', f'{self.style}.xsl')
    
    ####################### Define Values and attributes #######################

    def _init_values(self):
        """
        Method that defines the value objects for the Record.  This should
        call the super of this method, then use self._add_value to create new Value objects.
        Note that the order values are defined matters
        when build_model is called!!!
        """
        
        self._add_value('str', 'key',
                        valuerequired = True,
                        description = 'the UUID4 key for the record.')
        
        self._add_value('str', 'url',
                        modelpath = 'URL',
                        description = 'a URL where the record can be found.')
        
        self._add_value('str', 'method',
                        allowedvalues = (
                            'dynamic',
                            'static',
                            'box'),
                        description = 'the relaxation method used')
        
        self._add_value('str', 'standing',
                        allowedvalues = (
                            'good',
                            'bad'),
                        description = 'the current standing for the crystal with "good" records being unique untransformed entries')
        
        self._add_value('str', 'potential_LAMMPS_key',
                        valuerequired = True,
                        modelpath = 'potential-LAMMPS.key',
                        description = 'the UUID4 for the LAMMPS potential implementation used')
        
        self._add_value('str', 'potential_LAMMPS_id',
                        valuerequired = True,
                        modelpath = 'potential-LAMMPS.id',
                        description = 'the unique ID for the LAMMPS potential implementation used')
        
        self._add_value('str', 'potential_LAMMPS_url',
                        modelpath = 'potential-LAMMPS.URL',
                        description = 'a URL where the LAMMPS potential implementation record can be found')
        
        self._add_value('str', 'potential_key',
                        valuerequired = True,
                        modelpath = 'potential-LAMMPS.potential.key',
                        description = 'the UUID4 for the potential used')
        
        self._add_value('str', 'potential_id',
                        valuerequired = True,
                        modelpath = 'potential-LAMMPS.potential.id',
                        description = 'the unique ID for the potential used')
        
        self._add_value('str', 'potential_url',
                        modelpath = 'potential-LAMMPS.potential.URL',
                        description = 'a URL where the potential record can be found')
        
        self._add_value('float', 'temperature',
                        modelpath = 'phase-state.temperature',
                        metadatakey = 'T (K)',
                        unit = 'K',
                        description = 'the relaxation temperature')
        
        self._add_value('float', 'pressure_xx',
                        modelpath = 'phase-state.pressure-xx',
                        metadatakey = 'Pxx (GPa)',
                        unit = 'GPa',
                        description = 'tThe xx component of the relaxation pressure.')
        
        self._add_value('float', 'pressure_yy',
                        modelpath = 'phase-state.pressure-yy',
                        metadatakey = 'Pyy (GPa)',
                        unit = 'GPa',
                        description = 'the yy component of the relaxation pressure.')
        
        self._add_value('float', 'pressure_zz',
                        modelpath = 'phase-state.pressure-zz',
                        metadatakey = 'Pzz (GPa)',
                        unit = 'GPa',
                        description = 'the zz component of the relaxation pressure.')
        
        self._add_value('float', 'pressure_xy',
                        modelpath = 'phase-state.pressure-xy',
                        metadatakey = 'Pxy (GPa)',
                        unit = 'GPa',
                        description = 'the xy component of the relaxation pressure.')
        
        self._add_value('float', 'pressure_xz',
                        modelpath = 'phase-state.pressure-xz',
                        metadatakey = 'Pxz (GPa)',
                        unit = 'GPa',
                        description = 'the xz component of the relaxation pressure.')
        
        self._add_value('float', 'pressure_yz',
                        modelpath = 'phase-state.pressure-yz',
                        metadatakey = 'Pyz (GPa)',
                        unit = 'GPa',
                        description = 'the yz component of the relaxation pressure.')
        
        self._add_value('str', 'family',
                        valuerequired = True,
                        modelpath = 'system-info.family',
                        description = 'the ID of the reference structure family')
        
        self._add_value('str', 'family_url',
                        modelpath = 'system-info.family-URL',
                        description = 'a URL where the reference structure family can be found')
        
        self._add_value('str', 'parent_key',
                        valuerequired = True,
                        modelpath = 'system-info.parent_key',
                        description = 'the UUID4 key of the parent relaxation calculation record')
        
        self._add_value('str', 'parent_url',
                        modelpath = 'system-info.parent-URL',
                        description = 'a URL where the parent relaxation calculation record can be found')
        
        self._add_value('strlist', 'symbols',
                        valuerequired = True,
                        modelpath = 'system-info.symbol',
                        description = 'the potential element symbols in the crystal')
        
        self._add_value('str', 'composition',
                        valuerequired = True,
                        modelpath = 'system-info.composition',
                        description = 'the composition of the unit cell')
        
        self._add_value('str', 'crystalfamily',
                        modelpath = 'system-info.cell.crystal-family',
                        description = 'the crystal family of the unit cell: cubic, hexagonal, etc')
        
        self._add_value('int', 'natypes',
                        valuerequired = True,
                        modelpath = 'system-info.cell.natypes',
                        description = 'the number of unique atom types in the unit cell')
        
        self._add_value('float', 'a',
                        valuerequired = True,
                        modelpath = 'system-info.cell.a',
                        description = 'the a lattice constant for the unit cell')
                        
        self._add_value('float', 'b',
                        valuerequired = True,
                        modelpath = 'system-info.cell.b',
                        description = 'the b lattice constant for the unit cell')
                        
        self._add_value('float', 'c',
                        valuerequired = True,
                        modelpath = 'system-info.cell.c',
                        description = 'the c lattice constant for the unit cell')
                        
        self._add_value('float', 'alpha',
                        valuerequired = True,
                        modelpath = 'system-info.cell.alpha',
                        description = 'the alpha lattice angle for the unit cell')
                        
        self._add_value('float', 'beta',
                        valuerequired = True,
                        modelpath = 'system-info.cell.beta',
                        description = 'the beta lattice angle for the unit cell')
                        
        self._add_value('float', 'gamma',
                        valuerequired = True,
                        modelpath = 'system-info.cell.gamma',
                        description = 'the gamma lattice angle for the unit cell')
                        
        self._add_value('system_model', 'ucell',
                        valuerequired = True,
                        modelpath = "atomic-system",
                        description = 'the unit cell for the crystal')
                        
        self._add_value('float', 'potential_energy',
                        modelpath = 'potential-energy',
                        metadatakey = 'Epot (eV/atom)',
                        unit = 'eV',
                        description = 'the per-atom potential energy')
                        
        self._add_value('float', 'cohesive_energy',
                        modelpath = 'cohesive-energy',
                        metadatakey = 'Ecoh (eV/atom)',
                        unit = 'eV',
                        description = 'the per-atom cohesive energy')

    def set_ucell_attributes(self):
        """
        auto sets the symbols, composition, crystalfamily, natypes, a, b, c,
        alpha, beta, and gamma class attributes based on the current ucell.
        """
        self.symbols = self.ucell.symbols
        self.composition = self.ucell.composition
        self.crystalfamily = self.ucell.box.identifyfamily()
        self.natypes = self.ucell.natypes
        self.a = self.ucell.box.a
        self.b = self.ucell.box.b
        self.c = self.ucell.box.c
        self.alpha = self.ucell.box.alpha
        self.beta = self.ucell.box.beta
        self.gamma = self.ucell.box.gamma