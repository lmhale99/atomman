Updates
=======

Version 1.5.1
-------------

- **atomman.unitconvert** has been updated to use yabadaba.UnitConverter
  ensuring full compatibility between unit conversions with the two packages.
  Despite the code changing, atomman.unitconvert should still behave
  identically to before. **WARNING**: make sure yabadaba is at least version
  0.3.2 to fix a bug associated this change!!!

- **atomman.defect.GRIP** added that manages inputs for the grand-canonical
  interface predictor (GRIP) algorithm for grain boundary generation and
  relaxation.
  
- **atomman.defect.GrainBoundary.from_model** class method is added that reads
  in grain boundary configuration information from a grain_boundary record.

- **atomman.defect.GrainBoundary.dlat** method added that computes the lattice
  thickness for the grain boundary interface, which is used by GRIP.

- Miller vector support added to **atomman.library.record.GrainBoundary**.

- Some code cleanup and reorganization in **atomman.defect.Boundary.dlat**.

- **atomman.Box.d_hkl** method added for computing interplanar spacings.

- **atomman.mep.ISMPath** updated to use tqdm for the calculation status 
  progress rather than print statements.
  
- **atomman.defect.InterstitialSite** added, which uses Voronoi analysis to
  identify interstitial sites in a given system.

- **atomman.lammps.newseed** method added that generates a new
  LAMMPS-compatible random number seed.  **atomman.lammps.seed** method added
  that tests if a given int can be used as a LAMMPS random number seed.

- **atomman.lammps.run** tinkered with again related to trying to get the
  LAMMPS errors to be passed to python in an informative way.

- Minor doc and docstring updates for typos and syntax warnings.

Version 1.5.0
-------------

- Overhaul of the Record objects with the new yabadaba version that greatly
  reduces the code and better generalizes record interactions.

- Fixes to support numpy 2.0.

- **atomman.lammps.run** now has a partition parameter to support multi-replica
  simulations, like NEB and temperature-accelerated dynamics.  Updates (again)
  to hopefully print better error messages extracted from LAMMPS runs.

- **atomman.defect.SurfaceEnergyEstimator** and
  **atomman.defect.surface_energy_estimate()** added that provide quick
  estimates of high index free surface energies based on a few low index values.
  The former is a class for performing the calculations, and the latter a
  function.

- Bug fix for **ase_Atoms** dump to use scaled positions instead of absolute to
  ensure box and atoms remain aligned with the loss of box origin control.

- **atomman.tools.miller.tostring** method added for generating Miller indices
  strings for planes and vectors.  This serves as the inverse function of
  fromstring.

- **standardize_cell** dump style added that serves as a wrapper around the
  more generalized spglib.standardize_cell() method for identifying and
  normalizing systems to their corresponding unit cells.

- **Database.fetch_mp_crystal(s)** have been updated for the new Materials
  Project APIs.

Version 1.4.11
--------------

- The structure generation methods of **atomman.defect.Dislocation** now all
  have a center parameter that allows for the dislocation's position to be
  shifted.  This is useful for setting up NEB runs.

- **atomman.lammps.run** now supports lammps_command paths that contain spaces.

- **atomman.defect.Boundary** has been further improved and some supporting
  tools developed.  Not fully finalized yet.

- **neb_replica** dump style added that creates the atomic configurations used
  by LAMMPS for defining the final (and intermediate) replicas.
  
- More tools and operations related to Miller crystal vectors and planes have
  been added.  These should make it possible for future revisions of the
  defect generators to be simplified and made more uniform.

- Record classes have been updated to support a URL field for assigning
  persistent identification (PID) values once uploaded to a CDCS database.

Version 1.4.10
--------------

- **atomman.lammps.Log** now has more options during flatten allowing for only
  a subset of simulation runs to be merged together.  Bug fix related to type
  identification during flatten.

- **pdb** dump style for protein database file format.  Mostly useful for the
  plotting tools that natively interpret this format.

- Transformations added to support hexagonal conventional cells of trigonal
  systems. These 't1' and 't2' options are now available in
  **atomman.tools.miller.vector_conventional_to_primitive**, 
  **atomman.tools.miller.vector_primitive_to_conventional**,
  and the **conventional_to_primitive** and **primitive_to_conventional** dump
  styles.

- **atomman.tools.miller.fromstring** method added for interpreting formatted
  string representations of Miller vectors and planes as numpy arrays.

- From **atomman.defect**, **FreeSurface**, **StackingFault** and
  **Dislocation** now have fromrecord() and fromdatabase() methods.  The
  fromrecord() methods read in parameters from an associated defect parameter
  set record, and fromdatabase() allows for a matching defect parameter set
  record to be fetched from a database. 

- **atomman.defect.Dislocation** updated to recognize Miller-Bravais 4 index
  vectors and planes.

- All Record classes now have a "database" attribute that sets a default
  database to associate with the record.  This makes it possible for the Record
  classes to have methods that fetch additional content from the database as
  needed.

- **atomman.defect.Boundary** class added for generating phase boundary and
  grain boundary systems.  The basic core of the generator is there and works
  for cubic systems.  To be improved in the next few versions as time allows.

- **atomman.plot** has new methods to allow for interactive 3D plots when
  working in a Jupyter environment.  Substantial changes to these plotting
  tools likely in subsequent versions of atomman.

Version 1.4.9
-------------

- **atomman.thermo.RDF** updates to support more general use.

- Handling of record queries updated to be consistent with yabadaba 0.2.0.

Version 1.4.8
-------------

- **atomman.thermo** module added that provides tools and reference models
  related to thermodynamic calculations.

- **atomman.tools.vect_angle** updated to allow for comparisons of multiple
  vector values at once.

- **conventional_to_primitive** and **primitive_to_conventional** dump styles
  added that convert between standard conventional and primitive unit cell
  settings.  These were created to make it easy to generate compatible cells
  in the two settings for performing crystal and lattice vector operations.

- **atomman.tools.miller.vector_conventional_to_primitive** and 
  **atomman.tools.miller.vector_primitive_to_conventional** transformation
  operations for face-centered basis settings have been changed to coincide
  with the new dump methods above.

- **atomman.defect.VolterraDislocation** classes now support specifying the m
  and n dislocation axes using str 'x', 'y', 'z' values.  The dislocation
  transformation basis code is now integrated into the class as well.

- **atomman.defect.Dislocation** A new init parameter ucell_setting allows for
  the lattice setting of the ucell to be given so that all lattice vectors can
  be used and explored rather than just integer crystal vectors relative to
  ucell.  Also, a dipole method has been added that allows for the generation
  of stable dislocation dipole atomic configurations.

- **atomman.defect.Strain** bug fix related to numpy changing how irregular
  arrays are represented.

- **atomman.plot.interpolate_contour** figsize option fixed and improved.

Version 1.4.7
-------------

- **atomman.load** Import of load styles is now fully modular. Note that as a 
  result of this, the individual **load_{style}** function calls have been
  removed.

- **atomman.load** and **atomman.dump** Fix to make the individual load
  and dump styles optional dependent on any additional package requirements.
  Versions 1.4.4-1.4.6 accidentally required that these optional packages be
  installed.

- **atomman.lammps.Log.flatten** Fix for a rare case associated with stopped
  LAMMPS simulations and restarts that occasionally resulted in columns
  inadvertently being interpreted as str values instead of float values. 

Version 1.4.6
-------------

- Import of dump styles is now fully modular and delayed.  This allows for
  new styles to be introduced that have additional package requirements without
  breaking all of atomman. The delayed loading also makes it possible for dump
  methods to call other dump or load methods without import errors. Note that 
  as a result of this, the individual **dump_{style}** function calls have been
  removed.

- **atomman.dump.primitive_cell** now works properly because of the above.

- Internal use of **atomman.tools.crystalsystem** functions changed to use the
  corresponding **atomman.Box** methods introduced in 1.4.4. instead.

- XSL and XSD files added for the defined record styles to better support
  integration of the reference records with potentials.nist.gov.


Version 1.4.5
-------------

- **atomman.defect.FreeSurface** now has a unique_shifts() method that uses
  crystal symmetry to filter out most symmetrically equivalent termination
  planes.

- **atomman.dump.primitive_cell** has been added that uses spglib to take an
  atomic system and return a new system corresponding to the identified
  primitive unit cell.

- **atomman.dump.phonopy_Atoms** is updated for newer versions of phonopy.

- **strain** methods have been added to **atomman.defect.VolterraDislocation**
  and its subclasses to provide the strain associated with the dislocation
  solutions.  Stress and displacement methods in **atomman.defect.Stroh** have
  been adjusted to improve calculation speed.

- **atomman.defect.pn_arctan_disldensity** added and 
  **atomman.defect.pn_arctan_disregistry** updated for consistency and new
  options.  These give classic Peierls-Nabarro dislocation width models.

- **load_table** method has been added to **atomman.cluster.BondAngleMap**
  allowing for the data generated by save_table to be read back in.

Version 1.4.4
-------------

- New methods added to **atomman.Box**

  - **reciprocal_vects** method added that computes the reciprocal lattice
    vectors associated with the Box's vectors.

  - **vector_crystal_to_cartesian** and **plane_crystal_to_cartesian** from
    **atomman.tools.miller** have been added as Box methods.

  - **position_relative_to_cartesian** replaces **atomman.System.unscale** and
    **position_cartesian_to_relative** replaces **atomman.System.scale** as the
    new operations are better named and only relate to Box information.

  - **identifyfamily**, **iscubic**, **ishexagonal**, **istetragonal**,
    **isrhombohedral**, **isorthorhombic**, **ismonoclinic** and
    **ismonoclinic** methods added from **atomman.tools**.

- **atomman.region.Plane** has new methods **operate**, **__eq__** and
  **isclose** for transforming and comparing Planes.

- **atomman.plot.interpolate_contour** reworked to allow 
  matplotlib.pyplot.axes to be passed through allowing for the color contour
  plots to be added on top of existing plots.  Options also added to turn off
  features.

- **atomman.defect.DifferentialDisplacement** has new plot_with_nye method that
  overlays Nye tensor color contours with the differential displacement plots.

- Overhaul of **atomman.library** operations reflecting that underlying
  database handling is now branched off into the separate yabadaba package.

- **Dislocation**, **FreeSurface**, **PointDefect** and **StackingFault**
  Record classes related to defect parameter sets have been moved from iprPy
  to **atomman.library.record**  This is to support future updates where these
  parameter sets can be directly passed to the defect generator classes.

- Typing hints added to all of atomman's code.

Version 1.4.3
-------------

- **atomman.library.Database** query options better ordered and default values
  updated.  retrieve methods added to allow for database records to be copied
  to local files.

- Bug fix for composition queries of relaxed and reference crystal records.

- Updates for KIM model handling due to updates with the potentials package.


Version 1.4.2
-------------

- **atomman.dump.pymatgen_Structure** updated for new pymatgen versions.

- **atomman.defect.DifferentialDisplacement** bug fix related to handling
  the atomcolor and atomcmap parameters.

- **atomman.tools** now imports aslist, iaslist, screen_input, uber_open_rmode,
  and atomic_info from potentials to remove duplicate code. 

- **atomman.library** various updates related to keeping record handling
  consistent with updates in potentials version 0.3.1.

Version 1.4.1
-------------

- **atomman.lammps.Log** bug fix for properly reading performance data
  for restart runs.

Version 1.4.0
-------------

- **atomman.library** and **atomman.settings** modules updated to reflect
  the reworked potentials package version 0.3.0.

- **atomman.load_lammps_potential** and **atomman.load** options 'prototype'
  and 'crystal' updated for the new library module.  load style
  'dft_reference' added.

- **atomman.lammps.Potential** now is a function that returns either a 
  potentials.record.PotentialLAMMPS or potentials.record.PotentialLAMMPSKIM
  object.

- **atomman.lammps.run** now has options for passing string input scripts
  rather than reading from files, and for turning off log file output.
  **atomman.lammps.checkversion** simplified due to the changes to run.

- **atomman.cluster.BondAngleMap** added for characterizing the three-body
  interactions as predicted by interatomic potentials.  

Version 1.3.7
-------------

- **atomman.dump.atom_data** bug fix for kim model potentials (now they work).

- **atomman.lammps.Log** now captures performance output.  A Simulation class
  is added to better represent each run/simulation.  The flatten method is 
  updated to return a new Simulation rather than overwriting the current data.
  New 'all' style added to flatten that will merge all runs without filtering
  out duplicate timesteps. 

- **atomman.defect.differential_displacement** option added to pass an existing
  matplotlib axes object to plot on rather than generating a new figure.  This
  allows for subplots to be constructed.

- **atomman.defect.DifferentialDisplacement** option added to pass an existing
  matplotlib axes object to plot on rather than generating a new figure.  This
  allows for subplots to be constructed.

- **atomman.mep** subpackage added for performing minimum energy pathway
  calculations. The contained Path classes represent an energy path and have
  built-in iteration methods.  The ISMPath uses the improved string method.

  **atomman.defect.GammaSurface** updated with path and build_path methods
  that help build mep Path objects for the GammaSurface.

  **atomman.defect.Strain** class added that improves upon the nye_tensor
  function.  The new class uses Cython for roughly a 2X speedup and is
  designed to be easier to use.  

  **atomman.defect.SDVPN** The sign of tau used by stress_energy with
  fullstress=False is flipped to correspond to the behavior of 
  stress_energy with fullstress=True.  New parameter added allowing for
  additional kwargs to be passed to the underlying scipy.optimize.minimize(). 

Version 1.3.6
-------------

- **atomman.tools.atomic_info** updated for recently assigned element names
  and to be more lenient for isotopes.

- **atomman.dump.atom_data** updated to support using kim commands for kim
  model potentials.

- **atomman.dump.lammps_commands** added - NOT DEBUGGED FOR 
  NON-CUBIC/ORTHORHOMBIC SYSTEMS!

Version 1.3.5
-------------

- **atomman.defect.GammaSurface** updates and fixes related to the units
  parameters for the plotting methods.

- **atomman.defect.SDVPN** bug fixes related to model() generation, loading,
  and the units parameters for the plotting methods.

- **atomman.Settings** is now a renaming/import of potentials.Settings. 

Version 1.3.4
-------------

- **atomman.defect.Dislocation** class added that handles the generation of
  dislocation monopole and periodic array of dislocation atomic configurations
  in a more user-friendly interface than the previous functions.

- **atomman.region.PlaneSet** class added that allows for a region/shape to be
  defined using a list of planes.  This allows for the construction of 
  multi-faceted and/or open-ended shapes.

- **atomman.Box.planes** changed so that the order of the planes returned is
  consistent with the underlying indices.

- **atomman.build_lammps_potential** inherited from potentials package.

Version 1.3.3
-------------

- **atomman.Settings** class added that inherits from the corresponding class
  in the potentials package.  This makes it possible for atomman to access the
  same local directory of records as the potentials package.
  
- **atomman.library** module added that extends the corresponding module from
  the potentials package to include support for crystal_prototype and 
  relaxed_crystal records.

- **atomman.load_lammps_potential** added that loads LAMMPS potential
  information and downloads parameter files from the NIST Interatomic
  Potentials Repository.

- **atomman.load_prototype** and **atomman.load_crystal** load options added
  that allow for new Systems to be generated based on crystal_prototype and
  relaxed_crystal records in the NIST Interatomic Potentials Repository.

- **atomman.defect.GammaSurface** class updated so that the RBF interpolated
  energies are smoothed across the periodic cell boundaries.

- Fix to keep the code compatible with Python 3.6 (which broke in version
  1.3.2)

Version 1.3.2
-------------

- **System.r0** added which finds the shortest interatomic spacing.

- **System.rotate** made more robust.

- **atomman.tools.miller.plane_crystal_to_cartesian** added that identifies
  the Cartesian normal associated with a crystallographic plane.

- **atomman.lammps.Potential** made consistent with
  potentials.LAMMPSPotential.  Upcoming versions of atomman will have
  potentials as a requirement eliminating the duplication: (this class will
  simply be a renaming of the class from potentials).

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
  dislocation solutions based on the slip plane, line direction and Burgers
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
  more plotting options than the old differential_displacement function while
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