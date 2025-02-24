# coding: utf-8

# Standard Python imports
from typing import Optional, Tuple

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record

class CrystalPrototype(Record):
    """
    Class for representing crystal_prototype records that describe common
    crystal prototypes.
    """
    
    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'crystal_prototype'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'crystal-prototype'

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
        self._add_value('str', 'commonname', valuerequired=True, 
                        modelpath='name')
        self._add_value('str', 'prototype', valuerequired=True)
        self._add_value('str', 'pearson', valuerequired=True,
                        modelpath='Pearson-symbol')
        self._add_value('str', 'strukturbericht', valuerequired=True,
                        modelpath='Strukturbericht')
        self._add_value('int', 'sg_number', valuerequired=True,
                        modelpath='space-group.number')
        self._add_value('str', 'sg_hm', valuerequired=True,
                        modelpath='space-group.Hermann-Maguin')
        self._add_value('str', 'sg_schoenflies', valuerequired=True,
                        modelpath='space-group.Schoenflies')
        self._add_value('base', 'wykoff', valuerequired=True,
                        modelpath='space-group.Wykoff', metadatakey=False)
        self._add_value('str', 'crystalfamily', valuerequired=True,
                        modelpath='system-info.cell.crystal-family')
        self._add_value('int', 'natypes', valuerequired=True,
                        modelpath='system-info.cell.natypes')
        self._add_value('system_model', 'ucell', valuerequired=True,
                        modelpath="atomic-system")

    @property
    def defaultname(self) -> Optional[str]:
        """str: The name to default to, usually based on other properties"""
        return self.id
