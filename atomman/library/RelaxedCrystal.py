from copy import deepcopy
from DataModelDict import DataModelDict as DM

import atomman.unitconvert as uc
from .. import System
from ..tools import identifyfamily
class RelaxedCrystal():
    
    def __init__(self, model):
        self.load(model=model)

    def __str__(self):
        return f'Relaxed Crystal {self.key}'

    @property
    def key(self):
        return self.__key

    @property
    def method(self):
        return self.__method
    
    @property
    def standing(self):
        return self.__standing
    
    @property
    def family(self):
        return self.__family

    @property
    def parent_key(self):
        return self.__parent_key

    @property
    def potential_LAMMPS_id(self):
        return self.__potential_LAMMPS_id

    @property
    def potential_LAMMPS_key(self):
        return self.__potential_LAMMPS_key

    @property
    def potential_id(self):
        return self.__potential_id

    @property
    def potential_key(self):
        return self.__potential_key

    @property
    def cohesive_energy(self):
        return self.__cohesive_energy

    @property
    def ucell(self):
        return self.__ucell

    def load(self, model):
        
        self.__model = model = DM(model)
        
        crystal = model['relaxed-crystal']
        self.__key = crystal['key']
        self.__method = crystal['method']
        self.__standing = crystal['standing']

        self.__family = crystal['system-info']['family']
        self.__parent_key = crystal['system-info']['parent_key']
        
        self.__potential_LAMMPS_id = crystal['potential-LAMMPS']['id']
        self.__potential_LAMMPS_key = crystal['potential-LAMMPS']['key']
        self.__potential_id = crystal['potential-LAMMPS']['potential']['id']
        self.__potential_key = crystal['potential-LAMMPS']['potential']['key']
        self.__cohesive_energy = uc.value_unit(crystal['cohesive-energy'])
        
        self.__ucell = System(model=model)

    def asmodel(self):
        
        return deepcopy(self.__model)

    def asdict(self):
        """
        Converts the structured content to a simpler dictionary.
        
        Returns
        -------
        dict
            A dictionary representation of the record's content.
        """
        params = {}
        params['key'] = self.key
        params['method'] = self.method
        params['standing'] = self.standing
        params['family'] = self.family
        params['parent_key'] = self.parent_key
        
        params['potential_LAMMPS_id'] = self.potential_LAMMPS_id
        params['potential_LAMMPS_key'] = self.potential_LAMMPS_key
        params['potential_id'] = self.potential_id
        params['potential_key'] = self.potential_key
        
        params['cohesive_energy'] = self.cohesive_energy

        params['ucell'] = ucell = self.ucell
        params['crystal_family'] = identifyfamily(ucell.box)
        params['natypes'] = ucell.natypes
        params['symbols'] = ucell.symbols
        
        params['a'] = ucell.box.a
        params['b'] = ucell.box.b
        params['c'] = ucell.box.c
        params['alpha'] = ucell.box.alpha
        params['beta'] = ucell.box.beta
        params['gamma'] = ucell.box.gamma
        params['natoms'] = ucell.natoms


        
        return params