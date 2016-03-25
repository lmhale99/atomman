AtomMan
=======

Atomistic Manipulation Toolkit
 
AtomMan: the Atomistic Manipulation Toolkit is a Python library for 
creating, representing, manipulating, and analyzing large-scale atomic 
systems of atoms. The focus of the package is to facilitate the rapid design 
and development of simulations that are fully documented and easily adaptable 
to new potentials, configurations, etc.  

The code has no requirements that limit which systems it can be used on, i.e.
it should work on Linux, Mac and Windows computers.

The latest release can be installed using pip::

    pip install atomman

The code and all documentation is hosted on GitHub and can be directly 
downloaded at: `https://github.com/usnistgov/atomman`_.  

Documentation and examples in the form of Jupyter Notebooks for the various 
components of the code can be found by clicking on the links in the description 
below.

1. The core of atomman features:

    - Classes for representing an atomic system:
        
        + `Box`_ class defining a general parallelopid with instant conversions
          to/from lattice vector, unit cell, and LAMMPS-style parameter 
          representations.
          
        + `Atoms`_ class that stores the per-atom positions and properties for
          a collection of atoms. The underlying structure of the class uses 
          Numpy arrays and views to efficiently store the values and allow for 
          fast calculations.
          
        + `System`_ class that combines an Atoms object and a System object, 
          along with handling periodic boundaries.
        
    - Full `unit conversion`_ capabilities.

2. The lammps module provides tools for setting up, running, and analyzing 
   LAMMPS simulations directly from Python.  
   
    - `General LAMMPS`_ functionality, such as creating systems in LAMMPS,
      writing input files, reading/writing atom and dump files, calling LAMMPS
      to run, and recieving LAMMPS errors and log file thermo data.
     
    - `Potential`_ class for reading LAMMPS-potential data models.  This class
      and the data models allow for the modular treatment of any LAMMPS 
      implemented interatomic potentials. In other words, it allows for 
      potentials to be easily swapped.
      
    - `atom_data`_ load() and dump() functions allow for reading/writing of 
      LAMMPS atom coordination files. These functions work with any LAMMPS 
      atom_style and units options.
      
    - `atom_dump`_ load() and dump() functions allow for reading/writing of 
      LAMMPS dump coordination files. Additional JSON files are created that 
      define the proper conversion of the dump files to Python, and provide 
      metadata information to the terms (units, complete name, data structure, 
      etc).      

3. The defects module the following components:

    - `point`_ defect creation. Add vacancies, substitutionals, positional 
      interstitials, and dumbbell interstitials to a system.
      
    - `Stroh`_ method calculation. Solve the anisotropic elasticity solution
      for a perfectly straight dislocation of any character using the Stroh
      method. Once solved, the elastic energy coefficient and position-based
      displacements and stresses can be obtained. Dislocation monopole systems
      can be constructed using the Stroh displacements. (Soon...)
      
    - `dd`_: differential displacement plot for evaluating the core structure
      of dislocations. (Soon...)
    
    - `nye`_ tensor calculation for evaluating the core structure of 
      dislocations. (Soon...)

4. The models module has tools for reading in Systems:

    - `crystal`_ reads in a atomic-system data model file and returns a System. 
      (documentation soon...)
    
    - `cif`_ has two functions for reading CIF crystallographic files. 
      (documentation soon...)
    
        + cif_cell() reads a CIF file and returns a System and list of 
          elements.
          
        + cif_load() reads a CIF file and returns a DataModelDict 
          representation of all of the lines.
          
5. The convert module has tools for converting Systems to other codes:

    - `ase_Atoms`_: from_ase_Atoms() and to_ase_Atoms() convert between 
      ase.Atoms classes and atomman.System classes. NOTE: ase has requirements
      that are not compatible with Windows systems. (documentation soon...)
    
    - `pymatgen_Structure`_: from_pymatgen_Structure() and 
      to_pymatgen_Structure() convert between 
      pymatgen.Structure classes and atomman.System classes. NOTE: pymatgen has 
      requirements that are not compatible with Windows systems. (documentation soon...)
      
6. The tools module collects all of the additional features that don't have 
   homes anywhere else.
   
   - `ElasticConstants`_ class for interacting and transforming representations
     of the elastic constant tensor.
   
   - `nlist`_ function for constructing neighbor lists. (documentation soon...)
   
   - `plot tools`_ that help make pretty looking plots of properties. (soon ...)
   
   - `other tools`_ for system manipulation, calculations, etc. (documentation soon...)
   
.. _https://github.com/usnistgov/atomman: https://github.com/usnistgov/atomman
.. _Box: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.Box.ipynb
.. _Atoms: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.Atoms.ipynb
.. _System: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.System.ipynb
.. _unit conversion: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.unitconvert.ipynb
.. _General LAMMPS: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.lammps.ipynb
.. _Potential: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.lammps.Potential.ipynb
.. _atom_data: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.lammps.atom_data.ipynb
.. _atom_dump: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.lammps.atom_dump.ipynb
.. _point: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.defect.point.ipynb
.. _ElasticConstants: https://github.com/usnistgov/atomman/blob/master/Notebooks/atomman.tools.ElasticConstants.ipynb

