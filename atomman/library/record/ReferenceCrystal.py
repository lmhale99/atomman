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
        
        self._add_value('str', 'key',
                        valuerequired = True,
                        description = 'the UUID4 key for the record')

        self._add_value('str', 'id',
                        valuerequired = True,
                        description = 'the unique ID for the record')

        self._add_value('str', 'url',
                        modelpath = 'URL',
                        description = 'a URL where the record can be found')

        self._add_value('str', 'sourcename',
                        valuerequired = True,
                        modelpath = 'source.name',
                        description = 'the ID for the crystal used by the source site')
        
        self._add_value('str', 'sourcelink',
                        valuerequired = True,
                        modelpath = 'source.link',
                        description = 'a URL for the site where the crystal was retrieved')
        
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
