AtomMan
-------

Atomic Manipulation Toolkit 

This package is designed to allow easy manipulation of large-scale atomistic systems and the rapid development and use of property calculation metrics. 

Atomistic systems are represented using three basic classes: Atoms, Box, and System.  Atoms stores atomic properties using a numpy array to minimize memory requirements and allow for the properties to be used and set using numpy vectorization.  Box defines a generic parallelopid with methods for accessing different representations.  System combines an Atoms object and a Box object, along with the handling of periodic boundaries.

The lammps module has methods and classes to assist in the creation of LAMMPS 
input scripts and the exchange of information between LAMMPS and atomman.

The tools module contains numerous methods associated with analysis of the systems.

The model module contains tools associated with interacting with json/xml based data models.

The defect module contains routines for creating and analysing structural defects.

Demonstration Jupyter Notebooks for using the different classes, functions, modules and methods can be found in the contained Notebooks directory.