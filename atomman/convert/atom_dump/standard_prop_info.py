# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# atomman imports
from ...lammps import style

def standard_prop_info(lammps_units='metal'):
    
    lammps_unit = style.unit(lammps_units)
    
    return [{"prop_name": "a_id",
             "table_name": "id"},

            {"prop_name": "m_id",
             "table_name": "mol"},

            {"prop_name": "p_id",
             "table_name": "proc"},

            {"prop_name": "p_id_plus1",
             "table_name": "procp1"},

            {"prop_name": "atype",
             "table_name": "type"},

            {"prop_name": "element",
             "table_name": "element"},

            {"prop_name": "mass",
             "table_name": "mass",
             "unit": lammps_unit['mass']},

            {"prop_name": "pos",
             "table_name": ["x", "y", "z"],
             "unit": lammps_unit['length']},

            {"prop_name": "spos",
             "table_name": ["xs", "ys", "zs"],
             "unit": 'scaled'},

            {"prop_name": "upos",
             "table_name": ["xu", "yu", "zu"],
             "unit": lammps_unit['length']},

            {"prop_name": "supos",
             "table_name": ["xsu", "ysu", "zsu"],
             "unit": 'scaled'},

            {"prop_name": "boximage",
             "table_name": ["ix", "iy", "iz"],
             "unit": 'scaled'},

            {"prop_name": "velocity",
             "table_name": ["vx", "vy", "vz"],
             "unit": lammps_unit['velocity']},

            {"prop_name": "force",
             "table_name": ["fx", "fy", "fz"],
             "unit": lammps_unit['force']},

            {"prop_name": "charge",
             "table_name": "q",
             "unit": lammps_unit['charge']},

            {"prop_name": "mu",
             "table_name": ["mux", "muy", "muz"],
             "unit": lammps_unit['dipole']},

            {"prop_name": "mu_mag",
             "table_name": "mu",
             "unit": lammps_unit['dipole']},

            {"prop_name": "radius",
             "table_name": "radius",
             "unit": lammps_unit['length']},

            {"prop_name": "diameter",
             "table_name": "diameter",
             "unit": lammps_unit['length']},

            {"prop_name": "ang_velocity",
             "table_name": ["omegax", "omegay", "omegaz"],
             "unit": lammps_unit['ang-vel']},

            {"prop_name": "ang_momentum",
             "table_name": ["angmomx", "angmomy", "angmomz"],
             "unit": lammps_unit['ang-mom']},

            {"prop_name": "torque",
             "table_name": ["tqx", "tqy", "tqz"],
             "unit": lammps_unit['force'] + '*' + lammps_unit['length']}]