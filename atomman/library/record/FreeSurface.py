# coding: utf-8

# Standard Python imports
from typing import Tuple

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class FreeSurface(Record):
    """
    Class for representing free_surface records, which collect the parameters
    necessary for atomman to generate a particular free surface.
    """

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'free_surface'

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
        return 'free-surface'

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
        self._add_value('str', 'family', valuerequired=True,
                        modelpath='system-family')
        self._add_value('str', 'family_url',
                        modelpath='system-family-URL')
        self._add_value('miller', 'hkl', valuerequired=True,
                        modelpath='calculation-parameter.hkl',
                        bracket='()')
        self._add_value('int', 'shiftindex',
                        modelpath='calculation-parameter.shiftindex')
        self._add_value('str', 'cutboxvector', valuerequired=True,
                        modelpath='calculation-parameter.cutboxvector',
                        defaultvalue='c', allowedvalues=['a', 'b', 'c'])

    
    @property
    def hkl_str(self) -> str:
        """str: The string representation of hkl"""
        return self.get_value('hkl').str

    @property
    def parameters(self) -> dict:
        """dict : Defect parameters for atomman structure generator"""
        p = {}
        p['hkl'] = self.hkl_str
        if self.shiftindex is not None:
            p['shiftindex'] = str(self.shiftindex)
        p['cutboxvector'] = self.cutboxvector

        return p