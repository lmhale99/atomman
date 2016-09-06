AtomMan
=======

Atomistic Manipulation Toolkit
 
AtomMan: the Atomistic Manipulation Toolkit is a Python library for 
creating, representing, manipulating, and analyzing large-scale atomic 
systems of atoms. The focus of the package is to facilitate the rapid design 
and development of simulations that are fully documented and easily adaptable 
to new potentials, configurations, etc.  The code has no requirements that 
limit which systems it can be used on, i.e. it should work on Linux, Mac and 
Windows computers.

Features:

1. Allows for efficient and fast calculations on millions of atoms, each with many freely defined per-atom properties.
2. Create dislocation monopoles and evaluate them with differential displacement and Nye tensor plots.
3. Generate point defects.
4. Call LAMMPS directly from Python and instantly retrieve the resulting data or LAMMPS error statement.
5. Easily convert systems to/from the other Python atomic environments of ASE and PyMatGen.
6. Can create systems based on CIF crystal structure files, and LAMMPS atom and dump files.
7. Built-in unit conversions.

Installation
------------

The atomman package is designed for Python 2.7. It makes heavy use of numpy, so
it's easiest to download a Python environment like Anaconda. 

The latest release can be installed using pip::

    pip install atomman

This pip command should install atomman and any other required packages, but
occasionally a requirement may have to be installed separately. The list of required packages are given below.

Alternatively, all code and documentation can be downloaded from GitHub. 
    
    - The stable releases are available at `https://github.com/usnistgov/atomman`_.
    
    - The working development versions are at `https://github.com/lmhale99/atomman`_.
    
Documentation
-------------

Tutorials and full reference documentation can be found on GitHub in the form 
of Jupyter Notebooks. This provides explanations as well as examples of 
functioning code. They can also be downloaded and used interactively.

The links below are to the tutorials for the most recent stable release:

    - `Introduction to atomman`_
    
    - `Interacting with LAMMPS`_
    
    
Required packages
-----------------

This is a list of the required Python packages

    - `xmltodict`_ converts XML files to Python dictionaries. Used by 
      DataModelDict.
    
    - `DataModelDict`_ class allowing for easy transformations between 
      XML/JSON/Python representations of structured data models.
      
    - `numericalunits`_ helps with unit conversions.  
      
    - `numpy`_, `scipy`_, and `matplotlib`_ Python scientific tools
    
Optional packages
-----------------

This is a list of additional Python packages that can add functionality

    - `diffpy.Structure`_ - CIF reader. Required for loading systems from
      CIF files.
    
    - `ase`_ the Atomic Simulation Environment for interacting with small 
      systems and DFT calculations. Required for converting to/from ase.Atoms 
      objects.
    
    - `pymatgen`_ the Python Materials Genomics package used by the Materials
      Project for DFT calculations. Required for converting to/from 
      pymatgen.Structure objects.
       
.. _https://github.com/usnistgov/atomman: https://github.com/usnistgov/atomman
.. _https://github.com/lmhale99/atomman: https://github.com/lmhale99/atomman
.. _Introduction to atomman: https://github.com/usnistgov/atomman/blob/master/docs/tutorial/1%20Basics.ipynb
.. _Interacting with LAMMPS: https://github.com/usnistgov/atomman/blob/master/docs/tutorial/2%20LAMMPS%20Functionality.ipynb
.. _xmltodict: https://github.com/martinblech/xmltodict
.. _DataModelDict: https://github.com/usnistgov/DataModelDict
.. _numericalunits: https://pypi.python.org/pypi/numericalunits
.. _numpy: http://www.numpy.org/
.. _scipy: https://www.scipy.org/
.. _matplotlib: http://matplotlib.org/
.. _diffpy.Structure: http://www.diffpy.org/diffpy.Structure/
.. _ase: https://wiki.fysik.dtu.dk/ase/
.. _pymatgen: http://pymatgen.org/






