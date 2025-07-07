# coding: utf-8

# Standard Python imports
from typing import Tuple

import numpy as np

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class GrainBoundary(Record):
    """
    Class for representing grain_boundary records, which collect the parameters
    necessary for atomman to generate a particular grain boundary.
    """

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'grain_boundary'

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
        return 'grain-boundary'

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
        self._add_value('str', 'type')
        self._add_value('str', 'axis')
        self._add_value('float', 'misorientation-angle')
        self._add_value('miller', 'auvw1',
                        modelpath='calculation-parameter.auvw1')
        self._add_value('miller', 'buvw1',
                        modelpath='calculation-parameter.buvw1')
        self._add_value('miller', 'cuvw1',
                        modelpath='calculation-parameter.cuvw1')
        self._add_value('miller', 'auvw2',
                        modelpath='calculation-parameter.auvw2')
        self._add_value('miller', 'buvw2', 
                        modelpath='calculation-parameter.buvw2')
        self._add_value('miller', 'cuvw2',
                        modelpath='calculation-parameter.cuvw2')        
        self._add_value('str', 'cutboxvector', valuerequired=True,
                        modelpath='calculation-parameter.cutboxvector',
                        defaultvalue='c', allowedvalues=['a', 'b', 'c'])
        self._add_value('str', 'cellsetting', valuerequired=True,
                        modelpath='calculation-parameter.cellsetting',
                        defaultvalue='p',
                        allowedvalues=['p', 'i', 'f', 'a', 'b', 'c', 't', 't1', 't2'])
    
    @property
    def uvws1(self):
        """numpy.ndarray: All three rotation vectors for grain 1"""
        return np.array([self.auvw1, self.buvw1, self.cuvw1])

    @property
    def uvws2(self):
        """numpy.ndarray: All three rotation vectors for grain 2"""
        return np.array([self.auvw2, self.buvw2, self.cuvw2])

    @property
    def parameters(self) -> dict:
        """dict : Defect parameters for atomman structure generator"""
        p = {}
        p['uvws1'] = self.uvws1
        p['uvws2'] = self.uvws2
        p['cutboxvector'] = self.cutboxvector
        if self.cellsetting is not None:
            p['conventional_setting'] = self.cellsetting

        return p
    