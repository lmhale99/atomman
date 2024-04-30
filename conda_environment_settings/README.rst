Conda environment settings
--------------------------

This folder contains conda environment.yml files for creating/updating conda
environments for using atomman.

- **atomman.environment.yml** includes atomman and all optional dependencies
  and serves as a quick means of building a new environment for using atomman.

- **atomman_3_*.environment.yml** define environments with only the required
  dependencies for atomman (no optional packages) and are tied to a specific
  older Python version.  These are primarily used to build wheel distributions.



