# coding: utf-8
# Standard Python imports
from typing import Optional, Tuple

from scipy.spatial.transform import Rotation

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from DataModelDict import DataModelDict as DM

# Local imports
from .Boundary import Boundary
from .TiltGrainBoundaryHelper import TiltGrainBoundaryHelper
from .. import System, Box
from ..tools import approx_rational
from ..library.record.GrainBoundary import GrainBoundary as GBRecord


class GrainBoundary(Boundary):
    """
    Class for generating systems to investigate grain boundaries. In comparison
    to the more generic Boundary class, GrainBoundary only accepts one unit cell,
    it limits configurations to zero-strain states, and it includes extra methods
    and attributes specific to grain boundary configurations.
    """

    def __init__(self,
                 ucell: System,
                 uvws1,
                 uvws2,
                 conventional_setting: str = 'p',
                 cutboxvector='c',
                 maxmult: int = 10):
        """
        Initialize a grain boundary configuration.

        Parameters
        ----------
        ucell : atomman.System
            The reference unit cell to use for both grains.
        uvws1 : array-like object
            The three Miller(-Bravais) crystal vectors of ucell to use for
            orienting the first grain such that each crystal vector will be
            aligned with one of the box vectors of the final configuration.
        uvws2 : array-like object
            The three Miller(-Bravais) crystal vectors of ucell to use for
            orienting the second grain such that each crystal vector will be
            aligned with one of the box vectors of the final configuration.
        conventional_setting : str, optional
            Specifies which conventional lattice setting that ucell is in.
            The default value of 'p' takes ucell to be primitive, in which case
            the uvws1 and uvws2 values must be integers.  This must be specified
            in order to access non-integer lattice vectors for the uvws. 
        cutboxvector : str, optional
            Indicates which of the three box vectors that the boundary will
            be placed along. Default value is 'c'.
        maxmult : int, optional
            The max integer multiplier to use for the zero strain check.
        """
        # Call Boundary's init using the single ucell
        super().__init__(ucell, ucell, uvws1, uvws2, 
                         conventional_setting1=conventional_setting,
                         conventional_setting2=conventional_setting,
                         cutboxvector=cutboxvector,
                         zerostrain=True, maxmult=maxmult)
        
        # Find the transformation for rcell1 to rcell2
        rvects1 = self.transform_p1_to_r1.apply(self.ucell.box.vects)
        rvects2 = self.transform_p2_to_r2.apply(self.ucell.box.vects)
        self.transform = np.linalg.lstsq(rvects1, rvects2, rcond=None)[0].T
        self.__transform_r1_to_r2 = Rotation.from_matrix(np.linalg.lstsq(rvects1, rvects2, rcond=None)[0].T)

    @classmethod
    def from_model(cls,
                   ucell: System,
                   record: GBRecord,
                   maxmult: int = 10):
        """
        Allows for the grain boundary settings to be loaded from a grain_boundary
        record rather than manually specifying the inputs.

        Parameters
        ----------
        ucell : atomman.System
            The reference unit cell to use for both grains.
        record : str, path, DataModelDict.DataModelDict or atomman.library.record.GrainBoundary
            A grain_boundary-style record containing the input settings to load.
            Can be given as a file path, str file contents, DataModelDict contents,
            or a GrainBoundary record object.
        maxmult : int, optional
            The max integer multiplier to use for the zero strain check.
        
        Returns
        -------
        atomman.defect.GrainBoundary
        """
        # Extract inputs from the grain_boundary record
        if not isinstance(record, GBRecord):
            record = GBRecord(model = record)
        obj = cls(ucell, maxmult=maxmult, **record.parameters)

        return obj

    @classmethod
    def symmetric_tilt(cls,
                       ucell: System,
                       axis_uvw: npt.ArrayLike,
                       plane1_hkl: Optional[npt.ArrayLike] = None,
                       in1_uvw: Optional[npt.ArrayLike] = None,
                       conventional_setting: str = 'p',
                       ref_hkl: Optional[npt.ArrayLike] = None,
                       ref_uvw: Optional[npt.ArrayLike] = None,
                       cutboxvector: str = 'c',
                       maxindex: int = 20,
                       maxmult: int = 10,
                       tol: float = 1e-8):
        """
        Initialize for symmetric tilt grain boundary configurations.

        Parameters
        ----------
        ucell : atomman.System
            The unit cell to use for both top and bottom grains.
        axis_uvw : array-like
            The Miller [uvw] or Miller-Bravais [uvtw] crystal vector of ucell
            that serves as the tilt axis.  The tilt axis is a common vector in
            the grain boundary plane that is shared by both grains.
        in1_uvw : array-like object
            The first grain's in-plane vector given as a Miller or
            Miller-Bravais crystal vector.
        ref_uvw : array-like object or None, optional
            A reference in-plane vector that when crossed with the tilt axis
            identifies the grain boundary plane associated with the 0 degree
            misorientation configuration.  If given, in2 is found such
            that the two grain boundary planes have the same angle wrt the
            reference plane.  If None (default), then in2 is selected based on
            producing the smallest misorientation angle, which usually but not
            always corresponds to a symmetric tilt boundary.
        cutboxvector : str, optional
            Sets the alignment of the uvws to the final configuration's box
            vectors by specifying which of the box vectors is out-of-plane
            with respect to the grain boundary plane. Default value is 'c',
            which aligns the tilt axis with a, and in1 (and in2) with b.
        maxindex : int, optional
            The maximum absolute vector index to search over for the u,v,w
            values: u, v, and w can all independently vary from
            -maxindex to maxindex.  Default value is 20.
        maxmult : int, optional
            The max integer multiplier to use for the zero strain check if
            zerostrain is True.
        tol : float, optional
            A floating point tolerance used in finding the compatible primitive
            unit cell.  Available as a parameter in case issues are encountered.
            Default value is 1e-8.
        
        Returns
        -------
        atomman.defect.GrainBoundary
        """
        # Initialize helper
        helper = TiltGrainBoundaryHelper(ucell, axis_uvw,
                                         conventional_setting=conventional_setting,
                                         ref_uvw=ref_uvw, ref_hkl=ref_hkl, tol=tol)

        # Identify rotation uvws
        uvws1, uvws2 = helper.symmetric_uvws(plane1_hkl=plane1_hkl, in1_uvw=in1_uvw,
                                            cutboxvector=cutboxvector,
                                            maxindex=maxindex)
        obj = cls(ucell, uvws1, uvws2, conventional_setting=conventional_setting,
                  cutboxvector=cutboxvector)

        return obj

    @property
    def ucell(self) -> System:
        """atomman.System: The conventional reference unit cell for both grains. Alias for ucell1 and ucell2"""
        return self.ucell1

    @property
    def ucell_prim(self) -> System:
        """atomman.System: The primitive reference unit cell for both grains. Alias for ucell_prim1 and ucell_prim2"""
        return self.ucell_prim1
    
    @property
    def transform_r1_to_r2(self) -> Rotation:
        """scipy.spatial.transform.Rotation: The Cartesian rotation associated with rcell1 to rcell2"""
        return self.__transform_r1_to_r2
    
    @property
    def misorientation(self) -> float:
        """float : The misorientation angle for the grain boundary"""
        return np.linalg.norm(self.transform_r1_to_r2.as_rotvec(degrees=True))
    
    def sigma(self,
              tol: float = 1e-6) -> int:
        """Compute the sigma factor for the grain boundary"""

        t1 = self.transform_r1_to_r2.as_matrix()
        t2 = self.transform_r1_to_r2.inv().as_matrix()

        denominators = approx_rational(np.vstack([t1, t2]).flatten(), tol=tol)[1]

        # Sigma is the least common multiplier of the denominators
        return np.lcm.reduce(denominators)
    
    def dlat(self,
             precision: int = 6) -> float:
        """
        Computes the lattice thickness for the grain boundary by building a
        lattice system for the primitive cell and finding the distance between
        atoms perpendicular to the grain boundary plane.
        
        Parameters
        ----------
        precision : int, optional
            The rounding precision to use in identifying unique atom
            coordinates.  Default value is 6.

        Returns
        -------
        float
            The distance between two lattice atom positions in the direction
            perpendicular to the grain boundary plane.
        """
        # Construct a lattice system from the primitive unit cell and 1 atom
        latticesystem = System(box=Box(vects=self.ucell_prim.box.vects))
        
        # Rotate and supersize
        testsystem = latticesystem.rotate(self.uvws_prim1).supersize(2, 2, 2)
        
        # Find the difference between the two nearest atom positions perpendicular to the cut
        unique_z = sorted(list(set(testsystem.atoms.pos[:, self.cutindex].round(precision))))
        return unique_z[1] - unique_z[0]
    