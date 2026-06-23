# coding: utf-8
# Standard Python libraries
import io
from typing import Optional, Tuple, Union
import uuid

# https://github.com/usnistgov/DataModelDict
from DataModelDict import DataModelDict as DM

# https://github.com/usnistgov/yabadaba
from yabadaba.record import Record
from yabadaba import load_value

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

from scipy.interpolate import CubicSpline

__all__ = ['PartialStructureFactor']

class PartialStructureFactor(Record):

    def __init__(self,
                 model: Union[str, io.IOBase, DM, None] = None,
                 name: Optional[str] = None,
                 **kwargs):
        """
        Initializes a PartialStructureFactor Record object.
        
        Parameters
        ----------
        model : str, file-like object or DataModelDict, optional
            A JSON/XML data model for the content.
        name : str, optional
            The unique name to assign to the record.  If model is a file
            path, then the default record name is the file name without
            extension.
        system : str, optional
            The materials system that the partial structure factor was obtained
            from.  This is typically the list of all element symbols in the
            binary system and optionally any extra identifying information.
        elements : list, optional
            The two element symbols associated with the atomic pair interactions
            that the partial structure factor is for.
        q : array-like, optional
            The inverse radial coordinates.
        S : array-like, optional
            The partial structure factor, S(q), values for the pair of elements
            in the binary system.
        source : str, optional
            Information on the record's source, i.e. citation and/or creator
            details.
        """
        super().__init__(model=model, name=name, **kwargs)

    ########################## Basic metadata fields ##########################

    @property
    def style(self) -> str:
        """str: The record style"""
        return 'partial_structure_factor'

    @property
    def modelroot(self) -> str:
        """str: The root element of the content"""
        return 'partial-structure-factor'

    @property
    def xsl_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsl html transformer"""
        return ('dbliquid.xsl', 'partial-structure-factor.xsl')

    @property
    def xsd_filename(self) -> Tuple[str, str]:
        """tuple: The module path and file name of the record's xsd schema"""
        return ('dbliquid.xsd', 'partial-structure-factor.xsd')


    ############################# Define Values  ##############################

    def _init_values(self):
        """
        Method that defines the value objects for the Record.  This should
        call the super of this method, then use self._add_value to create new Value objects.
        Note that the order values are defined matters
        when build_model is called!!!
        """

        self._add_value('str', 'key', defaultvalue=str(uuid.uuid4()))
        self._add_value('str', 'id')
        self._add_value('longstr', 'source')
        self._add_value('longstr', 'system')
        self._add_value('strlist', 'elements')
        self._add_value('floatarray', 'q', valuerequired=True,
                        modelpath='partial-structure-factor-plot.q', unit='angstrom^-1')
        self._add_value('floatarray', 'S', valuerequired=True,
                        modelpath='partial-structure-factor-plot.S(q)')

    ########################## Additional settings ############################

    @property
    def defaultname(self) -> Optional[str]:
        """str: The name to default to, usually based on other properties"""
        return self.key

    ############################ Additional Methods ###########################

    def spline_S(self,
                 npoints: Optional[int] = None,
                 include_zero: bool = False,
                 inplace: bool = False):
        """
        Perform a cubic spline interpolation fit to the S(q) function.  One use
        for this is to interpolate missing intermediate points for S(q) tables
        that have irregularly spaced q values.

        Parameters
        ----------
        npoints : int or None, optional
            The number of q interpolation points to use in the range of q
            values.  If None (default) then the q values will be selected with
            a spacing equal to the smallest current q spacing.
        include_zero : bool, optional
            If True, then S(q) for q values less than given will be
            extrapolated by assuming S(q) smoothly goes to 0 at q=0.
        inplace: bool, optional
            If False (default), a new PartialStructureFactor object will be
            returned with the splined values.  If True, the values of the
            current object will be updated instead.

        Returns
        -------
        PartialStructureFactor
            If inplace is False, a new PartialStructureFactor object containing the
            interpolated values will be returned.
        """

        if include_zero:
            # Add S(q=0) = 0
            q = np.array([0.0] + self.q.tolist())
            S = np.array([0.0] + self.S.tolist())
            bc_type = ['clamped', 'not-a-knot']

        else:
            q = self.q
            S = self.S
            bc_type = ['not-a-knot', 'not-a-knot']

        # Fit so S =0 and dS/dq = 0 at q = 0
        spline = CubicSpline(q, S, bc_type=bc_type)

        if npoints is not None:
            if include_zero:
                newq = np.linspace(0.0, q[-1], npoints+1)[1:]
            else:
                newq = np.linspace(q[0], q[-1], npoints)

        else:
            Δq = np.min(q[1:] - q[:-1])
            if include_zero:
                newq = np.arange(Δq, q[-1] + Δq, Δq)
            else:
                newq = np.arange(q[0], q[-1] + Δq, Δq)

        newS = spline(newq)

        if inplace:
            self.q = newq
            self.S = newS

        else:
            return PartialStructureFactor(q=newq, S=newS, system=self.system,
                                          elements=self.elements, source=self.source)
    