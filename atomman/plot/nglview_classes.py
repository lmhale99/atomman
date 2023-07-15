# coding: utf-8
# Standard Python libraries
import uuid

try:
    from nglview.base_adaptor import Structure, Trajectory
    from nglview.adaptor import register_backend
    from nglview.widget import NGLWidget
except ModuleNotFoundError:
    has_nglview = False
else:
    has_nglview = True

if has_nglview:
    class AtommanStructure(Structure):
        def __init__(self, atomman_system, ext='pdb', params={}):
            super().__init__()
            self.path = ''
            self.ext = ext
            self.params = params
            self._atomman_system = atomman_system

        def get_structure_string(self):
            return self._atomman_system.dump('pdb')
            
    @register_backend('atomman')
    class AtommanTrajectory(Trajectory, Structure):
        """atommantraj adaptor."""

        def __init__(self, atomman_systems):
            """
            Initializes from a list of atomman.System objects
            """
            self.ext = 'pdb'
            self.params = {}
            self.trajectory = atomman_systems
            self.id = str(uuid.uuid4())

        def get_coordinates(self, index):
            return self.trajectory[index].pos

        def get_structure_string(self):
            return self.trajectory[0].dump('pdb')
            
        @property
        def n_frames(self):
            return len(self.trajectory)

    def show_atomman(atomman_system, **kwargs):
        """
        Examples
        --------
        >>> import nglview as nv
        >>> import atomman as am
        >>> ucell = am.load('prototype', 'A1--Cu--fcc', a=4.05, symbols='Al')
        >>> w = nv.show_atomman(ucell)
        """
        return NGLWidget(AtommanStructure(atomman_system), **kwargs)

    def show_atommantraj(atomman_systems, **kwargs):
        """Show trajectory using list of atomman.System objects.

        Examples
        --------
        >>> import nglview as nv
        >>> import atomman as am
        >>> systems = []
        >>> traj = Trajectory(nv.datafiles.ASE_Traj)
        >>> view = nv.show_asetraj(traj)
        >>> view.add_spacefill()
        >>> view # doctest: +SKIP
        """
        return NGLWidget(AtommanTrajectory(atomman_systems), **kwargs)

# Define dummy functions if no nglview installed
else:
    def show_atomman(atomman_system, **kwargs):
        raise ModuleNotFoundError('nglview is required by this operation')

    def show_atommantraj(atomman_systems, **kwargs):
        raise ModuleNotFoundError('nglview is required by this operation')