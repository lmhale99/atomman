# coding: utf-8

# Standard Python imports
from typing import Tuple

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class StackingFault(Record):
    """
    Record that collects input parameters for creating atomic configurations
    to evaluate a particular generalized stacking fault map.
    """

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'stacking_fault'

    @property
    def xsd_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsd schema"""
        return ('atomman.library.xsd', f'{self.style}.xsd')

    @property
    def xsl_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsl transformer"""
        return ('atomman.library.xsl', f'{self.style}.xsl')

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'stacking-fault'

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

        self._add_value('str', 'family',
                        valuerequired = True,
                        modelpath = 'system-family',
                        description = 'the ID of the reference structure family')
        
        self._add_value('str', 'family_url', 
                        modelpath = 'system-family-URL',
                        description = 'a URL where the family reference system record can be found')
        
        self._add_value('miller', 'hkl',
                        valuerequired = True,
                        modelpath = 'calculation-parameter.hkl',
                        bracket = '()',
                        description = 'the Miller(-Bravais) plane for the stacking fault')
        
        self._add_value('miller', 'a1vect_uvw',
                        valuerequired = True,
                        modelpath = 'calculation-parameter.a1vect_uvw',
                        bracket = '[]',
                        description = 'one of the two in-plane lattice vectors that defines the periodic 2D generalized fault cell')
        
        self._add_value('miller', 'a2vect_uvw',
                        valuerequired = True,
                        modelpath = 'calculation-parameter.a2vect_uvw',
                        bracket = '[]',
                        description = 'one of the two in-plane lattice vectors that defines the periodic 2D generalized fault cell')
        
        self._add_value('int', 'shiftindex',
                        modelpath = 'calculation-parameter.shiftindex',
                        description = 'the shift index parameter for identifying which atomic planes to slice between (i.e. this determines the termination plane)')
        
        self._add_value('str', 'cutboxvector',
                        valuerequired = True,
                        modelpath = 'calculation-parameter.cutboxvector',
                        defaultvalue = 'c',
                        allowedvalues = (
                            'a',
                            'b',
                            'c'),
                        description = 'indicates which system box vector is not in the stacking fault plane')

    
    @property
    def hkl_str(self) -> str:
        """str: The string representation of hkl"""
        return self.get_value('hkl').str

    @property
    def a1vect_uvw_str(self) -> str:
        """str: The string representation of a1vect_uvw"""
        return self.get_value('a1vect_uvw').str

    @property
    def a2vect_uvw_str(self) -> str:
        """str: The string representation of a2vect_uvw"""
        return self.get_value('a2vect_uvw').str

    @property
    def parameters(self) -> dict:
        """dict : Defect parameters for atomman structure generator"""
        p = {}
        p['hkl'] = self.hkl_str
        p['a1vect_uvw'] = self.a1vect_uvw_str
        p['a2vect_uvw'] = self.a2vect_uvw_str
        if self.shiftindex is not None:
            p['shiftindex'] = str(self.shiftindex)
        p['cutboxvector'] = self.cutboxvector

        return p
