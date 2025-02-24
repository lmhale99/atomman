# coding: utf-8

# Standard Python imports
from typing import Tuple

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class Dislocation(Record):
    """
    Class for representing dislocation records, which collect the parameters
    necessary for atomman to generate a particular dislocation type.
    """

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'dislocation'

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
        return 'dislocation'

    ####################### Define Values and attributes #######################

    def _init_values(self):
        """
        Method that defines the value objects for the Record.  This should
        call the super of this method, then use self._add_value to create new Value objects.
        Note that the order values are defined matters
        when build_model is called!!!
        """
        
        self._add_value('str', 'key', valuerequired=True)
        self._add_value('str', 'id',  valuerequired=True)
        self._add_value('str', 'url', modelpath='URL')
        self._add_value('str', 'family', valuerequired=True,
                        modelpath='system-family')
        self._add_value('str', 'family_url',
                        modelpath='system-family-URL')
        self._add_value('str', 'character', valuerequired=True,
                        allowedvalues=['screw', 'edge', 'mixed'])
        self._add_value('miller', 'slip_hkl', valuerequired=True,
                        modelpath='calculation-parameter.slip_hkl',
                        bracket='()')
        self._add_value('miller', 'line_uvw', valuerequired=True,
                        modelpath='calculation-parameter.line_uvw',
                        bracket='[]')
        self._add_value('miller', 'burgers_uvw', valuerequired=True,
                        modelpath='calculation-parameter.burgers_uvw',
                        bracket='[]')
        self._add_value('unitvector', 'm',
                        modelpath='calculation-parameter.m')
        self._add_value('unitvector', 'n',
                        modelpath='calculation-parameter.n')
        self._add_value('vector', 'shift',
                        modelpath='calculation-parameter.shift')
        self._add_value('bool', 'shiftscale',
                        modelpath='calculation-parameter.shiftscale')
        self._add_value('int', 'shiftindex',
                        modelpath='calculation-parameter.shiftindex')
   
    @property
    def slip_hkl_str(self) -> str:
        """str: The string representation of slip_hkl"""
        return self.get_value('slip_hkl').str
    
    @property
    def line_uvw_str(self) -> str:
        """str: The string representation of line_uvw"""
        return self.get_value('line_uvw').str

    @property
    def burgers_uvw_str(self) -> str:
        """str: The string representation of burgers_uvw"""
        return self.get_value('burgers_uvw').str

    @property
    def parameters(self) -> dict:
        """dict : Defect parameters for atomman structure generator"""
        p = {}
        p['slip_hkl'] = self.slip_hkl_str
        p['line_uvw'] = self.line_uvw_str
        p['burgers_uvw'] = self.burgers_uvw_str
        p['m'] = f"{self.m[0]} {self.m[1]} {self.m[2]}"
        p['n'] = f"{self.n[0]} {self.n[1]} {self.n[2]}"
        if self.shift is not None:
            p['shift'] = f"{self.shift[0]} {self.shift[1]} {self.shift[2]}"
        if self.shiftscale is not None:
            p['shiftscale'] = str(self.shiftscale)
        if self.shiftindex is not None:
            p['shiftindex'] = str(self.shiftindex)

        return p
