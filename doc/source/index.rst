=======================================
atomman: Atomistic Manipulation Toolkit
=======================================

Description
===========

AtomMan: the Atomistic Manipulation Toolkit is a Python library for
creating, representing, manipulating, and analyzing large-scale atomic
systems of atoms. The focus of the package is to facilitate the rapid design
and development of simulations that are fully documented and easily adaptable
to new potentials, configurations, etc.  The code has no requirements that
limit which systems it can be used on, i.e. it should work on Linux, Mac and
Windows computers.

Features:

1. Allows for efficient and fast calculations on millions of atoms, each with
   many freely defined per-atom properties.

2. Built-in tools for generating and analyzing crystalline defects, such as
   point defects, stacking faults, and dislocations.

4. Call LAMMPS directly from Python and instantly retrieve the resulting data
   or LAMMPS error statement.

5. Easily convert systems to/from the other Python atomic representations, such
   as ase.Atoms and pymatgen.Structure.

6. Can read and dump crystal structure information from a number of formats,
   such as LAMMPS data and dump files, and POSCAR.

7. Built-in unit conversions.

Installation
============

The atomman package is Python 3.7+ compatible.

The latest release can be installed using pip::

    pip install atomman

or, alternatively using conda and conda-forge::

    conda install atomman -c conda-forge

For Windows users, it is recommended to use an Anaconda distribution and use
conda to install numpy, scipy, matplotlib, pandas and cython prior to
installing atomman.

Alternatively, all code and documentation can be downloaded from GitHub.

- The stable releases are available at
  `https://github.com/usnistgov/atomman <https://github.com/usnistgov/atomman>`__.

- The working development versions are at
  `https://github.com/lmhale99/atomman <https://github.com/lmhale99/atomman>`__.

Tutorials
=========

.. toctree::
   :maxdepth: 2
   
   tutorial/index
    
    
Code Documentation
==================

.. toctree::
   :maxdepth: 3
   
   atomman

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`

.. _https://github.com/usnistgov/atomman: https://github.com/usnistgov/atomman
.. _https://github.com/lmhale99/atomman: https://github.com/lmhale99/atomman