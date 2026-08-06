# coding: utf-8

# Standard Python imports
from typing import Tuple

import numpy as np

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class GrainBoundary(Record):
    """
    Record that collects input parameters for creating atomic configurations
    of grain boundaries with a particular orientation.
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
                        description = 'a URL where the reference structure family can be found')
        
        self._add_value('str', 'type',
                        description = 'the category type of the grain boundary')

        self._add_value('str', 'axis',
                        description = 'the tilt/twist axis if there is one')

        self._add_value('float', 'misorientation',
                        modelpath = 'misorientation-angle',
                        description = 'the misorientation angle between the two grains')
        
        self._add_value('miller', 'auvw1',
                        modelpath = 'calculation-parameter.auvw1',
                        description = "the Miller-(Bravais) crystal vector of grain 1 that is aligned with the system's a box vector")
        
        self._add_value('miller', 'buvw1',
                        modelpath = 'calculation-parameter.buvw1',
                        description = "the Miller-(Bravais) crystal vector of grain 1 that is aligned with the system's b box vector")
        
        self._add_value('miller', 'cuvw1',
                        modelpath = 'calculation-parameter.cuvw1',
                        description = "the Miller-(Bravais) crystal vector of grain 1 that is aligned with the system's c box vector")
        
        self._add_value('miller', 'auvw2',
                        modelpath = 'calculation-parameter.auvw2',
                        description = "the Miller-(Bravais) crystal vector of grain 2 that is aligned with the system's a box vector")
        
        self._add_value('miller', 'buvw2', 
                        modelpath = 'calculation-parameter.buvw2',
                        description = "the Miller-(Bravais) crystal vector of grain 2 that is aligned with the system's b box vector")
        
        self._add_value('miller', 'cuvw2',
                        modelpath = 'calculation-parameter.cuvw2',
                        description = "the Miller-(Bravais) crystal vector of grain 2 that is aligned with the system's c box vector")
        
        self._add_value('str', 'cutboxvector',
                        valuerequired = True,
                        modelpath = 'calculation-parameter.cutboxvector',
                        defaultvalue = 'c',
                        allowedvalues = (
                            'a',
                            'b',
                            'c'),
                        description = 'indicates which system box vector is not in the grain boundary interface plane')
        
        self._add_value('str', 'cellsetting',
                        valuerequired = True,
                        modelpath = 'calculation-parameter.cellsetting',
                        defaultvalue = 'p',
                        allowedvalues = (
                            'p',
                            'i',
                            'f',
                            'a',
                            'b',
                            'c',
                            't',
                            't1',
                            't2'),
                        description = 'the crystal cell setting associated with the reference structure family (primitive, face-centered, etc) which allows for non-integer lattice vectors to be discovered')
    
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
    