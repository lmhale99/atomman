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
        
        self._add_value('str', 'key', valuerequired=True)
        self._add_value('str', 'url', modelpath='URL')
        self._add_value('str', 'method')
        self._add_value('str', 'standing')
        self._add_value('str', 'potential_LAMMPS_key',
                        valuerequired=True,
                        modelpath='potential-LAMMPS.key')
        self._add_value('str', 'potential_LAMMPS_id',
                        valuerequired=True,
                        modelpath='potential-LAMMPS.id')
        self._add_value('str', 'potential_LAMMPS_url',
                        modelpath='potential-LAMMPS.URL')
        self._add_value('str', 'potential_key', valuerequired=True,
                        modelpath='potential-LAMMPS.potential.key')
        self._add_value('str', 'potential_id', valuerequired=True,
                        modelpath='potential-LAMMPS.potential.id')
        self._add_value('str', 'potential_url',
                        modelpath='potential-LAMMPS.potential.URL')
        self._add_value('float', 'temperature',
                        #defaultvalue=0.0,
                        modelpath='phase-state.temperature',
                        metadatakey='T (K)', unit='K')
        self._add_value('float', 'pressure_xx',
                        #defaultvalue=0.0,
                        modelpath='phase-state.pressure-xx',
                        metadatakey='Pxx (GPa)', unit='GPa')
        self._add_value('float', 'pressure_yy',
                        #defaultvalue=0.0,
                        modelpath='phase-state.pressure-yy',
                        metadatakey='Pyy (GPa)', unit='GPa')
        self._add_value('float', 'pressure_zz',
                        #defaultvalue=0.0,
                        modelpath='phase-state.pressure-zz',
                        metadatakey='Pzz (GPa)', unit='GPa')
        self._add_value('float', 'pressure_xy',
                        #defaultvalue=0.0,
                        modelpath='phase-state.pressure-xy',
                        metadatakey='Pxy (GPa)', unit='GPa')
        self._add_value('float', 'pressure_xz',
                        #defaultvalue=0.0,
                        modelpath='phase-state.pressure-xz',
                        metadatakey='Pxz (GPa)', unit='GPa')
        self._add_value('float', 'pressure_yz',
                        #defaultvalue=0.0,
                        modelpath='phase-state.pressure-yz',
                        metadatakey='Pyz (GPa)', unit='GPa')
        self._add_value('str', 'family', valuerequired=True,
                        modelpath='system-info.family')
        self._add_value('str', 'family_url',
                        modelpath='system-info.family-URL')
        self._add_value('str', 'parent_key', valuerequired=True,
                        modelpath='system-info.parent_key')
        self._add_value('str', 'parent_url',
                        modelpath='system-info.parent-URL')
        self._add_value('strlist', 'symbols', valuerequired=True,
                        modelpath='system-info.symbol')
        self._add_value('str', 'composition', valuerequired=True,
                        modelpath='system-info.composition')
        self._add_value('str', 'crystalfamily',
                        modelpath='system-info.cell.crystal-family')
        self._add_value('int', 'natypes', valuerequired=True,
                        modelpath='system-info.cell.natypes')
        self._add_value('float', 'a', valuerequired=True,
                        modelpath='system-info.cell.a')
        self._add_value('float', 'b', valuerequired=True,
                        modelpath='system-info.cell.b')
        self._add_value('float', 'c', valuerequired=True,
                        modelpath='system-info.cell.c')
        self._add_value('float', 'alpha', valuerequired=True,
                        modelpath='system-info.cell.alpha')
        self._add_value('float', 'beta', valuerequired=True,
                        modelpath='system-info.cell.beta')
        self._add_value('float', 'gamma', valuerequired=True,
                        modelpath='system-info.cell.gamma')
        self._add_value('system_model', 'ucell', valuerequired=True,
                        modelpath="atomic-system")
        self._add_value('float', 'potential_energy',
                        modelpath='potential-energy',
                        metadatakey='Epot (eV/atom)', unit='eV')
        self._add_value('float', 'cohesive_energy',
                        modelpath='cohesive-energy',
                        metadatakey='Ecoh (eV/atom)', unit='eV')

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