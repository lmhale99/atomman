# coding: utf-8

# Standard Python imports
from typing import Tuple

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class PointDefectParameters(Record):
    """
    Class for representing the parameters section of the PointDefect record
    below.  Each parameter set describes a single point defect generation
    operation, with some point defects requiring multiple such operations.
    """
    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'point_defect_parameters'
   
    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'calculation-parameter'
    
    ####################### Define Values and attributes #######################

    def _init_values(self):
        """
        Method that defines the value objects for the Record.  This should
        call the super of this method, then use self._add_value to create new Value objects.
        Note that the order values are defined matters
        when build_model is called!!!
        """
        
        self._add_value('str', 'ptd_type', allowedvalues=['v', 'i', 's', 'db'])
        self._add_value('int', 'atype')
        self._add_value('vector', 'pos')
        self._add_value('vector', 'db_vect')
        self._add_value('bool', 'scale')

    @property
    def parameters(self) -> dict:
        """dict : Defect parameters for atomman structure generator"""
        p = {}
        p['ptd_type'] = self.ptd_type
        if self.atype is not None:
            p['atype'] = str(self.atype)
        if self.pos is not None:
            p['pos'] = f"{self.pos[0]} {self.pos[1]} {self.pos[2]}"
        if self.db_vect is not None:
            p['db_vect'] = f"{self.db_vect[0]} {self.db_vect[1]} {self.db_vect[2]}"
        if self.scale is not None:
            p['scale'] = str(self.scale)

        return p

class PointDefect(Record):
    """
    Class for representing point_defect records, which collect the parameters
    necessary for atomman to generate a particular point defect.
    """

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'point_defect'

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
        return 'point-defect'
    
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
        self._add_value('record', 'parameter_records', recordclass=PointDefectParameters,
                        modelpath='calculation-parameter',
                        metadatakey='parameters')
        
    @property
    def parameters(self) -> list:
        """list : Defect parameters for atomman structure generator"""
        allparameters = []
        for subset in self.parameter_records:
            allparameters.append(subset.parameters)
        
        return allparameters