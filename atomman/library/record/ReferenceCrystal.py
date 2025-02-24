# coding: utf-8

# Standard Python imports
from typing import Optional, Tuple

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class ReferenceCrystal(Record):
    """
    Class for representing reference_crystal records that provide the structure
    information for DFT relaxed crystal structures obtained from DFT databases.
    """

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'reference_crystal'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'reference-crystal'

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
        self._add_value('str', 'id', valuerequired=True)
        self._add_value('str', 'url', modelpath='URL')
        self._add_value('str', 'sourcename', valuerequired=True,
                        modelpath='source.name')
        self._add_value('str', 'sourcelink', valuerequired=True,
                        modelpath='source.link')
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

    @property
    def defaultname(self) -> Optional[str]:
        """str: The name to default to, usually based on other properties"""
        return self.id
       
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
