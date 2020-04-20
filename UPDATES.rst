Updates
=======

Version 1.3.2
-------------

- **System.r0** added which finds the shortest interatomic spacing.

- **System.rotate** made more robust.

- **atomman.tools.miller.plane_crystal_to_cartesian** added that identifies
  the Cartesian normal associated with a crystallographic plane.

- **atomman.lammps.Potential** made consistent with
  potentials.LAMMPSPotential.  Version 1.4 of atomman will have potentials
  as a requirement eliminating the duplication: (this class will simply be
  a renaming of the class from potentials).

- **atomman.lammps.LammpsError** error type added.

- **atomman.defect.dislocation_system_basis** and 
  **atomman.defect.dislocation_system_transform** functions added supporting
  the identification of dislocation system orientations based on
  material-specific parameters.  
  
- The "n" parameter in **atomman.defect.free_surface_basis** was renamed to
  maxindex consistency with the new dislocation_system functions.

- **atomman.defect.VolterraDislocation**, **atomman.defect.Stroh**,
  **atomman.defect.IsotropicVolterraDislocation**, and
  **atomman.defect.solve_volterra_dislocation** were updated by integrating in
  the dislocation_system functions. This makes it possible to now easily define
  dislocaiton solutions based on the slip plane, line direction and Burgers
  vector alone.
  
- **atomman.defect.dislocation_periodic_array** was updated to add an old_id
  parameter to the returned dislocation system making it easier to map the atoms
  in the defect system back to the perfect crystal base system used during
  construction.
  
- **atomman.defect.FreeSurface** class for generating free surface
  configurations from a unit cell and (hkl) plane was added.

- **atomman.defect.StackingFault** class completely rebuilt as a subclass of
  FreeSurface to make it easier to use, i.e. systems can be generated directly
  from unit cell, (hkl) and shift values.

- **atomman.defect.DifferentialDisplacement** class created. This class offers
  more ploting options than the old differential_displacement function while
  dividing the calculation and plotting into separate steps to make it easier
  to work with.

- **atomman.defect.SDVPN** class updated to allow for VolterraDislocation
  objects to be directly used as input parameters.  This makes it easier to
  work with as the transformations between dislocation orientations and gamma
  surface orientations can be automatically identified and handled.
  Additionally, solution summary and plotting tools incorporated into the
  class for convenience.

Version 1.3.1
-------------

- **Atoms.prop_atype** updated for new atype handling.

- **defect.GammaSurface** default plotting behavior improved.


Version 1.3.0
-------------

- **Support for Python < 3.6 removed.**  Python 2 support removed due to its
  imminent end at the new year.  Minimal version of 3.6 selected to take
  advantage of f-strings.

- **Atoms and System natype, atypes** behavior changed to allow for unassigned
  atype values and/or symbols.  Now, atype values must be > 0 and natypes =
  max(atype).  CAUTION: this could conceivably break backwards compatibility.

- **lammps.Potential** expanded.

    - **allsymbols** property added to support pair_styles that require all
      symbols to be listed in the pair_coeff lines even if they are not used.
    - **status** property added that indicates if the potential is known to
      have been superseded by a newer version or retracted for being invalid.
    - **pair_info** now supports an optional masses parameter for overriding
      default mass values.

- **load.atom_data** now recognizes image flags in the Atoms tables, and reads
  values from the Masses tables.  Parameter checking is performed allowing for
  more informative errors to be thrown.

- **dump.atom_data** updated to allow Potential objects to be passed directly,
  and for pair_info to be included in the generated info LAMMPS input lines.

- **System.masses** attribute added.  This is used for saving mass values from
  load.atom_data, and for overriding default Potential.masses values in
  dump.atom_data.

- **defect.dislocation_array** debugged, documented, and made consistent with
  Volterra solutions.

- **defect.IsotropicVolterraDislocation** displacements fixed and adjusted to
  predict displacements and stresses consistent with values from defect.Stroh.

- **defect.solve_volterra_dislocation** simplified to remove unnecessary 
  pre-check of elastic constants.

- **region** submodule added that allows for geometries in space to be defined
  and used to slice systems and per-atom properties.

- **Box** is now a subclass of region.Shape allowing it to be used for 
  region-based selection as well.

Version 1.2.8
-------------

- **defect.GammaSurface** support added for setting shift vectors using
  Miller-Bravais 4-term vectors.

- **tools.duplicates_allclose** added that identifies unique value sets
  based on absolute tolerances.

- **load('phonopy'), System.dump('phonopy')** bug fixes.

- **System.atoms_ix** compatibility checks changed and reduced from throwing
  an error to throwing a warning.

- **Atoms.extend and System.atoms_extend** methods added for adding atoms to
  existing Atoms/System objects.

Version 1.2.7
-------------

- **Atoms.model and Box.model** added to create/read data model 
  representations of the objects.

- **System.composition** added that returns string composition.

- **System.model, load('system_model'), System.dump('system_model')**
  data model format improved to capture all system information.

- **tools.Miller** functions for converting between Miller and Miller-Bravais
  crystal planes.

- **defect.GammaSurface** combining of multiple plots better supported.

- **defect.StackingFault** minimum r parameter added allowing all atoms to
  be at least a certain distance apart.

- **defect.free_surface_basis** added for identifying system orientations
  associated with free surface configurations.

Version 1.2.6
-------------

- **lammps.NEBLog** added for nudged elastic band calculation log files.

- **tools.Miller** transformations now all take float values and
  primitive-conventional cell conversions added.

- **Box.volume** bug fix to ensure returned volume is always positive.

- **defect.StackingFault** stacking fault configuration generator added.

- **nlist, dvect, dmag, defect.slip_vector** routines improved using Cython,
  alternate implementations of routines removed.

Version 1.2.5
-------------

- **Box.volume** parameter added.  Also, new class methods for initializing boxes
  based on crystal systems (cubic, hexagonal, etc.).

- **load('poscar')** now supports excess per-atom lines.

- **System.atoms_ix** added for indexing atoms at the system level.

- **defect.GammaSurface** reworked with improved design and features.

Version 1.2.4
-------------

- **Atoms.prop_atype()** added to allow properties to be assigned by prop_atype.

- **ElasticConstants.normalized_as()** and **ElasticConstants.is_normal()** added to
  force/check crystallographic symmetry of elastic tensors.

- **load('atom_data')** updated to support reading files containing # comments.

- **lammps.Potential** now supports specifying potentials with static charges.

- **defect.IsotropicVolterraDislocation** class added as **defect.Stroh** could not calculate
  isotropic solutions. Both classes are now children of **defect.VolterraDislocation**,
  and wrapper function **defect.solve_volterra_dislocation()** has been added.

- **defect.dislocation_array()** added that transforms a bulk system into a periodic array of
  dislocations, where the two system boundaries in the slip plane are periodic, and
  the third boundary is not.

- **defect.differential_displacement()** updated to provide users more options and control over
  the plots.

- MANIFEST.in corrected so non-code files should be properly copied during installation.

Version 1.2.3
-------------

- **load()** updated with more uniform parameters across the different styles.  
  Style 'phonopy_Atoms' added.

- **System.wrap()** made slightly more robust.

Version 1.2.2
-------------
- **System** scale/unscale bug fix.

- **defect.GammaSurface.model()** returned format improved for saving/loading results.

- **load('system_model')** updated with symbols parameter.

Version 1.2.1
-------------

- Corrections to setup.py for properly loading/building cython code.

Version 1.2.0
-------------

- Overhaul for Python 2/3 compatibility.

- Reorganization of code and renaming of some features.

- Cython routines added for dvect and neighbor list calculations.

- Improved documentation.