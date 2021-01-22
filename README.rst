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
------------

The atomman package is compatible with Python 3.6+.

The latest release can be installed using pip::

    pip install atomman

or using conda from the conda-forge channel::

    conda config --add channels conda-forge
    conda config --set channel_priority strict
    conda install <package-name>

For Windows users, it is recommended to use an Anaconda distribution and use
conda to install numpy, scipy, matplotlib, pandas and cython prior to
installing atomman.

Alternatively, all code and documentation can be downloaded from GitHub.

- The stable releases are available at
  `https://github.com/usnistgov/atomman <https://github.com/usnistgov/atomman>`__.

- The working development versions are at
  `https://github.com/lmhale99/atomman <https://github.com/lmhale99/atomman>`__.

Library setup
-------------
Starting with version 1.3.3 atomman uses the potentials Python package to
interact with the potentials.nist.gov database.  This allows for users to
easily search, discover and use the NIST-hosted interatomic potentials as well
as some computed properties associated with the potentials.  There are a few
setup options that can help you get the most out of this feature.

Note: the Python packages potentials, atomman and iprPy all share the same
settings file.  Updating the library settings in one package will carry over
to the other packages.

Changing the library directory
``````````````````````````````
Potential files and database records associated with the potentials can be
saved locally to a library directory.  The default location of the library
directory is <HOME>/.NISTPotentials/library. You can easily change the path
to a more accessible location if you wish::

    import atomman as am
    settings = am.Settings()
    settings.set_library_directory("NEWPATH")

Default behavior
````````````````
The atomman.Library class serves as the central means of interacting with the
records database.  By default, the class will search for matching records
first from the local library directory, then from the remote
potentials.nist.gov.  If records with the same name are found in both
locations, the local copy will be taken over the remote.  This makes it
possible for users to modify existing records and add their own user-defined
records. 

All of the Library methods that retrieve/load records and the load function
options that use the Library class have the following parameters

- **remote** (*bool*) indicates if the remote potentials.nist.gov will be
  searched.
- **local** (*bool*) indicates if the local library directory will be searched.
- **localpath** (*str*) allows for an alternate local library directory path
   to be searched.

Heavy-usage/offline behavior
````````````````````````````
For users who plan on extensive interactions with the database or who are
running on systems with limited internet availability, the following steps are
recommended:

1. Download/clone the github repository at
   https://github.com/lmhale99/potentials-library to the library directory.
   This github repository is a snapshot copy of all potential records at
   potentials.nist.gov and is kept (mostly) up to date.
2. In Python, load atomman.Settings() and set the default remote behavior to
   False::
 
    settings.set_remote(False)

   With this setting, the default value of remote parameter (above) for the
   Library class methods will be False. This eliminates the need for an
   internet connection (after step 1) and is typically much faster at
   retrieving records.
    
Documentation
-------------

Web-based documentation for the atomman package is available at
`https://www.ctcms.nist.gov/potentials/atomman <https://www.ctcms.nist.gov/potentials/atomman>`__.

Source code for the documentation can be found in the
`github doc directory <https://github.com/usnistgov/atomman/tree/master/doc/>`__.
The doc directory contains the information both as the source RestructuredText
files and as unformatted HTML. If you download a copy, you can view the HTML
version offline by

    cd {atomman_path}/doc/html
    python -m http.server

Then, opening localhost:8000 in a web browser.

The documentation consists of two main components:

1. **Tutorial Jupyter Notebooks:**
   `Online html version <https://www.ctcms.nist.gov/potentials/atomman/tutorial/index.html>`__,
   `Downloadable Notebook version <https://github.com/usnistgov/atomman/tree/master/doc/tutorial>`__.
   The tutorials starting with ##. provide a general overview/example of the
   various capabilities.  The tutorials starting with ##.#. give more detailed
   descriptions and list options available to the tools mentioned in the
   overview tutorials.

2. **Code Documentation:**
   `Online html version <https://www.ctcms.nist.gov/potentials/atomman/atomman.html>`__.
   This provides a rendering of the Python docstrings for the included
   functions and classes.


Optional packages
-----------------

This is a list of additional Python packages that are needed for some of the
optional features of the package.

- `diffpy.Structure <http://www.diffpy.org/diffpy.Structure/>`__:
  CIF reader. Required for loading systems from CIF files.

- `ase <https://wiki.fysik.dtu.dk/ase/>`__:
  The Atomic Simulation Environment for interacting with small systems
  and DFT calculations. Required for converting to/from ase.Atoms objects.

- `pymatgen <http://pymatgen.org/>`__:
  The Python Materials Genomics package used by the Materials
  Project for DFT calculations. Required for converting to/from
  pymatgen.Structure objects.

- `spglib <https://atztogo.github.io/spglib/python-spglib.html>`__:
  A Python interface to the spglib spacegroup analysis code.  spglib
  can be used to analyze and determine the spacegroup for an atomic system.
  Required for converting to/from spglib.cell objects.
