# coding: utf-8

# Standard Python imports
import io
from typing import Optional, Union, Tuple

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_query, load_value

# atomman imports
import atomman.unitconvert as uc
from ... import System

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

    def _init_value_objects(self) -> list:
        """
        Method that defines the value objects for the Record.  This should
        1. Call the method's super() to get default Value objects.
        2. Use yabadaba.load_value() to build Value objects that are set to
           private attributes of self.
        3. Append the list returned by the super() with the new Value objects.

        Returns
        -------
        value_objects: A list of all value objects.
        """
        value_objects = super()._init_value_objects()
        
        self.__key = load_value('str', 'key', self, valuerequired=True)
        self.__url = load_value('str', 'url', self, modelpath='URL')
        self.__method = load_value('str', 'method', self)
        self.__standing = load_value('str', 'standing', self)
        self.__potential_LAMMPS_key = load_value('str', 'potential_LAMMPS_key', self,
                                                 valuerequired=True,
                                                 modelpath='potential-LAMMPS.key')
        self.__potential_LAMMPS_id = load_value('str', 'potential_LAMMPS_id', self,
                                                valuerequired=True,
                                                modelpath='potential-LAMMPS.id')
        self.__potential_LAMMPS_url = load_value('str', 'potential_LAMMPS_url', self,
                                                modelpath='potential-LAMMPS.URL')
        self.__potential_key = load_value('str', 'potential_key', self, valuerequired=True,
                                         modelpath='potential-LAMMPS.potential.key')
        self.__potential_id = load_value('str', 'potential_id', self, valuerequired=True,
                                         modelpath='potential-LAMMPS.potential.id')
        self.__potential_url = load_value('str', 'potential_url', self,
                                           modelpath='potential-LAMMPS.potential.URL')
        self.__temperature = load_value('float', 'temperature', self,
                                        #defaultvalue=0.0,
                                        modelpath='phase-state.temperature',
                                        metadatakey='T (K)', unit='K')
        self.__pressure_xx = load_value('float', 'pressure_xx', self,
                                        #defaultvalue=0.0,
                                        modelpath='phase-state.pressure-xx',
                                        metadatakey='Pxx (GPa)', unit='GPa')
        self.__pressure_yy = load_value('float', 'pressure_yy', self,
                                        #defaultvalue=0.0,
                                        modelpath='phase-state.pressure-yy',
                                        metadatakey='Pyy (GPa)', unit='GPa')
        self.__pressure_zz = load_value('float', 'pressure_zz', self,
                                        #defaultvalue=0.0,
                                        modelpath='phase-state.pressure-zz',
                                        metadatakey='Pzz (GPa)', unit='GPa')
        self.__pressure_xy = load_value('float', 'pressure_xy', self,
                                        #defaultvalue=0.0,
                                        modelpath='phase-state.pressure-xy',
                                        metadatakey='Pxy (GPa)', unit='GPa')
        self.__pressure_xz = load_value('float', 'pressure_xz', self,
                                        #defaultvalue=0.0,
                                        modelpath='phase-state.pressure-xz',
                                        metadatakey='Pxz (GPa)', unit='GPa')
        self.__pressure_yz = load_value('float', 'pressure_yz', self,
                                        #defaultvalue=0.0,
                                        modelpath='phase-state.pressure-yz',
                                        metadatakey='Pyz (GPa)', unit='GPa')
        self.__family = load_value('str', 'family', self, valuerequired=True,
                                   modelpath='system-info.family')
        self.__family_url = load_value('str', 'family_url', self,
                                       modelpath='system-info.family-URL')
        self.__parent_key = load_value('str', 'parent_key', self, valuerequired=True,
                                       modelpath='system-info.parent_key')
        self.__parent_url = load_value('str', 'parent_url', self,
                                       modelpath='system-info.parent-URL')
        self.__symbols = load_value('strlist', 'symbols', self, valuerequired=True,
                                    modelpath='system-info.symbol')
        self.__composition = load_value('str', 'composition', self, valuerequired=True,
                                       modelpath='system-info.composition')
        self.__crystalfamily = load_value('str', 'crystalfamily', self,
                                          modelpath='system-info.cell.crystal-family')
        self.__natypes = load_value('int', 'natypes', self, valuerequired=True,
                                    modelpath='system-info.cell.natypes')
        self.__a = load_value('float', 'a', self, valuerequired=True,
                              modelpath='system-info.cell.a')
        self.__b = load_value('float', 'b', self, valuerequired=True,
                              modelpath='system-info.cell.b')
        self.__c = load_value('float', 'c', self, valuerequired=True,
                              modelpath='system-info.cell.c')
        self.__alpha = load_value('float', 'alpha', self, valuerequired=True,
                                  modelpath='system-info.cell.alpha')
        self.__beta = load_value('float', 'beta', self, valuerequired=True,
                                  modelpath='system-info.cell.beta')
        self.__gamma = load_value('float', 'gamma', self, valuerequired=True,
                                  modelpath='system-info.cell.gamma')
        self.__ucell = load_value('system_model', 'ucell', self, valuerequired=True,
                                  modelpath="atomic-system")
        self.__potential_energy = load_value('float', 'potential_energy', self,
                                             modelpath='potential-energy',
                                             metadatakey='Epot (eV/atom)', unit='eV')
        self.__cohesive_energy = load_value('float', 'cohesive_energy', self,
                                            modelpath='cohesive-energy',
                                            metadatakey='Ecoh (eV/atom)', unit='eV')

        value_objects.extend([
            self.__key, self.__url, self.__method, self.__standing,
            self.__potential_LAMMPS_key, self.__potential_LAMMPS_id, self.__potential_LAMMPS_url, 
            self.__potential_key, self.__potential_id,  self.__potential_url,
            self.__temperature, self.__pressure_xx, self.__pressure_yy, self.__pressure_zz,
            self.__pressure_xy,  self.__pressure_xz,  self.__pressure_yz,
            self.__family, self.__family_url, self.__parent_key, self.__parent_url,
            self.__symbols, self.__composition, self.__crystalfamily, self.__natypes, 
            self.__a, self.__b, self.__c, self.__alpha, self.__beta, self.__gamma,
            self.__ucell, self.__potential_energy, self.__cohesive_energy,
        ])

        return value_objects
   
    @property
    def key(self) -> str:
        """str : A UUID4 key assigned to the record"""
        return self.__key.value

    @key.setter
    def key(self, val: str):
        self.__key.value = val

    @property
    def url(self) -> Optional[str]:
        """str : A URL where a copy of the record can be found"""
        return self.__url.value

    @url.setter
    def url(self, val: str):
        self.__url.value = val

    @property
    def method(self) -> str:
        """str : Indicates the relaxation method used: box, static or dynamic"""
        return self.__method.value
    
    @method.setter
    def method(self, val: str):
        self.__method.value = val

    @property
    def standing(self) -> str:
        """str : 'good' or 'bad', with bad indicating it to be a duplicate or transformation"""
        return self.__standing.value
    
    @standing.setter
    def standing(self, val: str):
        self.__standing.value = val

    @property
    def potential_LAMMPS_key(self) -> str:
        """str : The key of the LAMMPS implementation used to relax the crystal"""
        return self.__potential_LAMMPS_key.value
    
    @potential_LAMMPS_key.setter
    def potential_LAMMPS_key(self, val: str):
        self.__potential_LAMMPS_key.value = val

    @property
    def potential_LAMMPS_id(self) -> str:
        """str : The id of the LAMMPS implementation used to relax the crystal"""
        return self.__potential_LAMMPS_id.value

    @potential_LAMMPS_id.setter
    def potential_LAMMPS_id(self, val: str):
        self.__potential_LAMMPS_id.value = val

    @property
    def potential_LAMMPS_url(self) -> Optional[str]:
        """str : A URL where a copy of the potential_LAMMPS record can be found"""
        return self.__potential_LAMMPS_url.value

    @potential_LAMMPS_url.setter
    def potential_LAMMPS_url(self, val: str):
        self.__potential_LAMMPS_url.value = val

    @property
    def potential_key(self) -> str:
        """str : The key of the potential model used to relax the crystal"""
        return self.__potential_key.value

    @potential_key.setter
    def potential_key(self, val: str):
        self.__potential_key.value = val

    @property
    def potential_id(self) -> str:
        """str : The id of the potential model used to relax the crystal"""
        return self.__potential_id.value

    @potential_id.setter
    def potential_id(self, val: str):
        self.__potential_id.value = val

    @property
    def potential_url(self) -> Optional[str]:
        """str : A URL where a copy of the potential model record can be found"""
        return self.__potential_url.value

    @potential_url.setter
    def potential_url(self, val: str):
        self.__potential_url.value = val

    @property
    def temperature(self) -> float:
        """float : The target temperature used during relaxation"""
        return self.__temperature.value

    @temperature.setter
    def temperature(self, val: float):
        self.__temperature.value = val

    @property
    def pressure_xx(self) -> float:
        """float : The target xx pressure component used during relaxation"""
        return self.__pressure_xx.value
    
    @pressure_xx.setter
    def pressure_xx(self, val: float):
        self.__pressure_xx.value = val

    @property
    def pressure_yy(self) -> float:
        """float : The target yy pressure component used during relaxation"""
        return self.__pressure_yy.value

    @pressure_yy.setter
    def pressure_yy(self, val: float):
        self.__pressure_yy.value = val

    @property
    def pressure_zz(self) -> float:
        """float : The target zz pressure component used during relaxation"""
        return self.__pressure_zz.value
    
    @pressure_zz.setter
    def pressure_zz(self, val: float):
        self.__pressure_zz.value = val

    @property
    def pressure_xy(self) -> float:
        """float : The target xy pressure component used during relaxation"""
        return self.__pressure_xy.value
    
    @pressure_xy.setter
    def pressure_xy(self, val: float):
        self.__pressure_xy.value = val

    @property
    def pressure_xz(self) -> float:
        """float : The target xz pressure component used during relaxation"""
        return self.__pressure_xz.value
    
    @pressure_xz.setter
    def pressure_xz(self, val: float):
        self.__pressure_xz.value = val

    @property
    def pressure_yz(self) -> float:
        """float : The target yz pressure component used during relaxation"""
        return self.__pressure_yz.value

    @pressure_yz.setter
    def pressure_yz(self, val: float):
        self.__pressure_yz.value = val

    @property
    def family(self) -> str:
        """str : The associated prototype/reference crystal id that the relaxed crystal is based on"""
        return self.__family.value

    @family.setter
    def family(self, val: str):
        self.__family.value = val

    @property
    def family_url(self) -> Optional[str]:
        """str : A URL where a copy of the family record can be found"""
        return self.__family_url.value

    @family_url.setter
    def family_url(self, val: str):
        self.__family_url.value = val

    @property
    def parent_key(self) -> str:
        """str : The key assigned to the record of the relaxation calculation used"""
        return self.__parent_key.value

    @parent_key.setter
    def parent_key(self, val: str):
        self.__parent_key.value = val

    @property
    def parent_url(self) -> Optional[str]:
        """str : A URL where a copy of the parent record can be found"""
        return self.__parent_url.value

    @parent_url.setter
    def parent_url(self, val: str):
        self.__parent_url.value = val

    @property
    def symbols(self) -> list:
        """list : The list of element model symbols"""
        return self.__symbols.value
    
    @symbols.setter
    def symbols(self, val: Union[str, list]):
        self.__symbols.value = val

    @property
    def composition(self) -> str:
        """str : The crystal's composition"""
        return self.__composition.value
    
    @composition.setter
    def composition(self, val: str):
        self.__composition.value = val

    @property
    def crystalfamily(self) -> str:
        """str : The crystal's system family"""
        return self.__crystalfamily.value
    
    @crystalfamily.setter
    def crystalfamily(self, val: str):
        self.__crystalfamily.value = val

    @property
    def natypes(self) -> int:
        """int : The number of atom types in the unit cell"""
        return self.__natypes.value
    
    @natypes.setter
    def natypes(self, val: int):
        self.__natypes.value = val

    @property
    def a(self) -> float:
        """float : The unit cell's a lattice parameter"""
        return self.__a.value
    
    @a.setter
    def a(self, val: float):
        self.__a.value = val

    @property
    def b(self) -> float:
        """float : The unit cell's b lattice parameter"""
        return self.__b.value
    
    @b.setter
    def b(self, val: float):
        self.__b.value = val

    @property
    def c(self) -> float:
        """float : The unit cell's c lattice parameter"""
        return self.__c.value
    
    @c.setter
    def c(self, val: float):
        self.__c.value = val

    @property
    def alpha(self) -> float:
        """float : The unit cell's alpha lattice angle"""
        return self.__alpha.value
    
    @alpha.setter
    def alpha(self, val: float):
        self.__alpha.value = val

    @property
    def beta(self) -> float:
        """float : The unit cell's beta lattice angle"""
        return self.__beta.value
    
    @beta.setter
    def beta(self, val: float):
        self.__beta.value = val

    @property
    def gamma(self) -> float:
        """float : The unit cell's gamma lattice angle"""
        return self.__gamma.value
    
    @gamma.setter
    def gamma(self, val: float):
        self.__gamma.value = val

    @property
    def ucell(self) -> System:
        """atomman.System : The unit cell system for the crystal"""
        return self.__ucell.value

    @ucell.setter
    def ucell(self, val: System):
        self.__ucell.value = val
    
    @property
    def potential_energy(self) -> float:
        """float : The measured per-atom potential energy"""
        return self.__potential_energy.value
    
    @potential_energy.setter
    def potential_energy(self, val: float):
        self.__potential_energy.value = val

    @property
    def cohesive_energy(self) -> float:
        """float : The computed per-atom cohesive energy"""
        return self.__cohesive_energy.value

    @cohesive_energy.setter
    def cohesive_energy(self, val: float):
        self.__cohesive_energy.value = val

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