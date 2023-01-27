# coding: utf-8
# cython: language_level=3

# Standard Python libraries
import warnings

# http://www.numpy.org/
import numpy as np

# atomman imports
from .. import NeighborList
from ..tools import axes_check, aslist

cimport cython

from libc.math cimport cos, pi, sqrt, fabs

class Strain():

    def __init__(self, system, neighbors=None, cutoff=None, p_vectors=None,
                 theta_max=27, axes=None, basesystem=None, baseneighbors=None):
        """
        Class initializer.  Allows for the current and reference state to be
        specified for performing the calculations.
        
        Parameters
        ----------
        system : atomman.System
            The atomic system to compute the per-atom strain properties and Nye
            tensor for.
        neighbors : atomman.NeighborList, optional
            The neighbor list associated with system to use.  Either neighbors
            or cutoff must be given, or system must have a neighbors attribute.
        cutoff : float
            Cutoff distance for computing a neighbor list for system.  Either
            neighbors or cutoff must be given, or system have a neighbors
            attribute.
        p_vectors : array-like object, optional
            List(s) of radial distance vectors between each atom and its nearest
            neighbors in a perfect crystal setting.  If one list of p_vectors is
            given, then it is applied to all atoms.
        axes : array-like object, optional
            3x3 array of right-handed orthogonal axes to transform the given p
            vectors.  Only needed if the orientation of the p vectors differs
            from the system's orientation.
        basesystem : atomman.System, optional
            A reference atomic system to use for constructing the p vectors. 
        baseneighbors : atomman.NeighborList, optional
            The neighbor list associated with basesystem to use. If basesystem
            is given, then either baseneighbors or cutoff must be given, or
            basesystem must have a neighbors attribute.
        theta_max : float, optional
            The maximum theta angle in degrees to use when searching for matches
            between p vectors and q vectors.  Optimum values are dependent on the
            crystal structure. Default value is 27, which is the original value
            used for fcc crystals.
        """
        
        self.__system = system
        
        # Neighbor list setup
        if neighbors is not None:
            assert cutoff is None, 'neighbors and cutoff cannot both be given'
            self.__neighbors = neighbors
        elif cutoff is not None:
            self.__neighbors = NeighborList(system=system, cutoff=cutoff)
        elif hasattr(system, 'neighbors'):
            self.__neighbors = system.neighbors
        else:
            raise ValueError('neighbors or cutoff is required')
        
        # p vector setup
        if basesystem is not None:
            assert p_vectors is None, 'basesystem and p_vectors cannot both be given'
            self.build_p_vectors(basesystem, neighbors=baseneighbors, cutoff=cutoff)
        elif p_vectors is not None:
            self.set_p_vectors(p_vectors, axes=axes)
        else:
            self.__p_vectors = None
    
        self.theta_max = theta_max
        self.clear_properties()
    
    @property
    def system(self):
        """atomman.System: The system the properties are being computed for."""
        return self.__system
    
    @property
    def p_vectors(self):
        """numpy.NDArray: The per-atom sets of ideal atom positions."""
        return self.__p_vectors
    
    @property
    def theta_max(self):
        """float: The maximum angle in degrees to include in the p-q vector pairings."""
        return self.__theta_max
    
    @theta_max.setter
    def theta_max(self, value):
        if value <= 180 and value > 0:
            self.__theta_max = value
    
    @property
    def neighbors(self):
        """atomman.NeighborList: The list of neighbors for system."""
        return self.__neighbors
    
    @property
    def G(self):
        """numpy.NDArray : The computed per-atom lattice correspondence tensor"""
        if self.__G is None:
            self.solve_G()
        return self.__G
    
    @property
    def strain(self):
        """numpy.NDArray : The computed per-atom strain tensor"""
        if self.__strain is None:
            self.__strain = strain_c(self.G)
        return self.__strain

    @property
    def invariant1(self):
        """numpy.NDArray : The computed per-atom first strain invariant"""
        if self.__invariant1 is None:
            self.__invariant1 = invariant1_c(self.strain)
        return self.__invariant1
        
    @property
    def invariant2(self):
        """numpy.NDArray : The computed per-atom second strain invariant"""
        if self.__invariant2 is None:
            self.__invariant2 = invariant2_c(self.strain)
        return self.__invariant2
        
    @property
    def invariant3(self):
        """numpy.NDArray : The computed per-atom third strain invariant"""
        if self.__invariant3 is None:
            self.__invariant3 = invariant3_c(self.strain)
        return self.__invariant3

    @property
    def rotation(self):
        """numpy.NDArray : The computed per-atom rotation tensor"""
        if self.__rotation is None:
            self.__rotation = rotation_c(self.G)
        return self.__rotation
        
    @property
    def angularvelocity(self):
        """numpy.NDArray : The computed per-atom angular velocity"""
        if self.__angularvelocity is None:
            self.__angularvelocity = angularvelocity_c(self.rotation)
        return self.__angularvelocity

    @property
    def nye(self):
        """numpy.NDArray : The computed per-atom Nye tensor"""
        if self.__nye is None:
            self.solve_nye()
        return self.__nye

    def clear_properties(self):
        """
        Clears all computed properties. Allows for the values to be recomputed
        using different settings.
        """
        self.__G = None
        self.__strain = None
        self.__invariant1 = None
        self.__invariant2 = None
        self.__invariant3 = None
        self.__angularvelocity = None
        self.__rotation = None
        self.__nye = None

    def save_to_system(self, properties=None):
        """
        Saves the computed per-atom strain properties to the given system.

        Parameters
        ----------
        properties : str or list, optional
            One or more properties.  If not given, will include strain,
            invariant1, invariant2, invariant3, angularvelocity and nye.
        """
        defaultkeys = ['strain', 'invariant1', 'invariant2', 'invariant3',
                       'angularvelocity', 'nye']
        allkeys = ['G', 'rotation'] + defaultkeys

        if properties is None:
            properties = defaultkeys
        else:
            properties = aslist(properties)
        for p in properties:
            assert p in allkeys, 'unknown property ' + p
            self.system.atoms.view[p] = getattr(self, p)

    def asdict(self, properties=None):
        """
        Returns a dictionary containing the computed per-atom strain properties.
        This corresponds to the results from the old nye_tensor function.

        Parameters
        ----------
        properties : str or list, optional
            One or more properties.  If not given, will include strain,
            invariant1, invariant2, invariant3, angularvelocity and nye.

        Returns
        -------
        dict
            Containing each of the computed properties.
        """
        defaultkeys = ['strain', 'invariant1', 'invariant2', 'invariant3',
                       'angularvelocity', 'nye']
        allkeys = ['G', 'rotation'] + defaultkeys
        results = {}
        if properties is None:
            properties = defaultkeys
        else:
            properties = aslist(properties)
        
        results = {}
        for p in properties:
            assert p in allkeys, 'unknown property ' + p
            results[p] = getattr(self, p)
        
        return results

    def set_p_vectors(self, p_vectors, axes=None):
        """
        Sets the p vectors for each atom according to a single p vector set to use
        for all atoms or a list of p vector sets defined for each atom individually.

        Parameters
        ----------
        p_vectors : array-like object
            List(s) of radial distance vectors between each atom and its nearest
            neighbors in a perfect crystal setting.  If one list of p_vectors is
            given, then it is applied to all atoms.
        axes : array-like object, optional
            3x3 array of right-handed orthogonal axes.  If given, will be used to
            transform the p_vectors before computing the Nye tensor.
        """
        system = self.system
        
        # Broadcast a single p_vectors list to all atoms
        if len(p_vectors) == 1:
            p_vectors = np.broadcast_to(p_vectors, (system.natoms, len(p_vectors[0]), 3))
        elif len(p_vectors) != system.natoms:
            p_vectors = np.broadcast_to(p_vectors, (system.natoms, len(p_vectors), 3))
        else:
            for i in range(len(p_vectors)):
                p_vectors[i] = np.asarray(p_vectors[i])
                if p_vectors[i].ndim == 1:
                    p_vectors[i] = np.array([p_vectors[i]])
            p_vectors = np.asarray(p_vectors)

        # Transform p_vectors if axes is given
        if axes is not None:
            p_vectors = np.inner(p_vectors, axes_check(axes))
            
        self.__p_vectors = p_vectors
        
    def build_p_vectors(self, basesystem, neighbors=None, cutoff=None):
        """
        Builds the p vectors for each atom based on a reference system.

        Parameters
        ----------
        basesystem : atomman.system
            The base/reference system to use.  This should be a defect-free
            perfect crystal system with atom ids directly corresponding to atoms
            in any system that you want to analyze with the Nye tensor.
        neighbors : atomman.NeighborList, optional
            The neighbor list associated with system to use.  Either neighbors
            or cutoff must be given, or system must have a neighbors attribute.
        cutoff : float
            Cutoff distance for computing a neighbor list for system.  Either
            neighbors or cutoff must be given, or system have a neighbors
            attribute.
        """

        # Neighbor list setup
        if neighbors is not None:
            assert cutoff is None, 'neighbors and cutoff cannot both be given'
        elif cutoff is not None:
            neighbors = NeighborList(system=basesystem, cutoff=cutoff)
        elif hasattr(basesystem, 'neighbors'):
            neighbors = basesystem.neighbors
        else:
            raise ValueError('neighbors or cutoff is required')

        p = []
        for i in range(basesystem.natoms):
            p.append(basesystem.dvect(i, neighbors[i]))

        self.__p_vectors = np.asarray(p, dtype=object)

    def solve_G(self, theta_max=None):
        """
        Computes the lattice correspondence tensor, G, by comparing all neighbor
        vectors q in the current system to ideal reference neighbor vectors p.
        
        Parameters
        ----------
        theta_max : float, optional
            Allows for the theta_max value set during the class initialization
            to be changed.  The maximum theta angle in degrees to use when
            searching for matches between p vectors and q vectors.  Optimum
            values are dependent on the crystal structure. 
        """
        # p vector setup
        p_vectors = self.p_vectors
        if p_vectors is None:
            raise ValueError('Cannot solve until p_vectors are set')

        # Theta_max setup
        if theta_max is not None:
            self.theta_max = theta_max
        cos_theta_max = cos(self.theta_max * pi / 180.0)

        # system and neighbors setup
        system = self.system
        neighbors = self.neighbors
        maxcoord = neighbors.coord.max()
        natoms = system.natoms

        # Initialize variables
        self.clear_properties()
        G = np.empty((natoms, 3, 3))
        P = np.empty((maxcoord, 3))
        Q = np.empty((maxcoord, 3))

        # Loop over all atoms
        for i in range(natoms):
            # Get actual p and q for atom i
            p = np.asarray(p_vectors[i], dtype=float)
            q = system.dvect(i, neighbors[i])
            
            # Get the matched P and Q
            n = match_pq(p, q, cos_theta_max, P, Q)
            
            # Compute lattice correspondence tensor, G, from P and Q
            if n == 0:
                G[i] = np.identity(3)
                warnings.warn('An atom lacks pair sets. Check neighbor list')
            else:
                G[i] = np.linalg.lstsq(Q[:n], P[:n], rcond=None)[0]

        self.__G = G

    def solve_nye(self):
        """
        Computes the Nye tensor based on the lattice correspondence tensor G
        and the current atomic positions.
        """
        G = self.G
        cdef double[:,:,::1] G_view = G
        system = self.system
        nlist = self.neighbors.nlist
        cdef long long[:,::1] nlist_view = nlist
        
        cdef Py_ssize_t i, x, y, z
        
        # Variable setup and shortcuts
        natoms = G.shape[0]
        cdef long long maxcoord = 0
        for i in range(natoms):
            if maxcoord < nlist_view[i, 0]:
                maxcoord = nlist_view[i, 0]
        
        nye = np.empty((natoms, 3, 3))
        cdef double [:,:,::1] nye_view = nye
        gradG = np.empty((3, 3, 3))
        cdef double [:,:,::1] gradG_view = gradG
        
        gG = np.empty((3,3))
        cdef double [:,::1] gG_view = gG
        
        dG = np.empty((maxcoord, 3, 3))
        cdef double [:,:,::1] dG_view = dG
        
        Q = np.empty((maxcoord, 3))
        cdef double [:,::1] Q_view = Q
        
        # Construct the gradient tensor of G, gradG for each atom
        for i in range(natoms):
            
            c = nlist_view[i, 0]
            Q[:c] = system.dvect(i, nlist[i, 1:c+1])
            
            dG_c(G_view, nlist_view, dG_view, i)
            
            for x in range(3):
                gG[:] = np.linalg.lstsq(Q[:c], dG[:c,x,:], rcond=None)[0]
                
                for y in range(3):
                    for z in range(3):
                        gradG_view[x,y,z] = gG_view[z,y]

            nye_c(gradG_view, nye_view, i)
            
        self.__nye = nye

@cython.boundscheck(False)
@cython.wraparound(False) 
cdef match_pq(double[:,::1] p, 
              double[:,::1] q,
              double cos_theta_max,
              double[:,::1] P, 
              double[:,::1] Q):
    """
    Takes p and q neighbor vectors and finds corresponding pairs P and Q based
    on the angles beween vectors (and interatomic spacing).
    """
    # Parameter initialization
    pnum = p.shape[0]
    qnum = q.shape[0]
    cdef long long[:] qp_pairs = np.empty(qnum, dtype='int64')
    cdef double [:] pmag = np.empty(pnum)
    cdef double [:] qmag = np.empty(qnum)
    
    cdef double r1 = 1e16
    cdef double cos_theta_min, jrad, krad
    
    cdef Py_ssize_t j, k, n, x, y, z
    
    # Compute magnitudes of q and set all qp_pairs to -1
    for j in range(qnum):
        qmag[j] = sqrt(q[j,0] * q[j,0] + q[j,1] * q[j,1] + q[j,2] * q[j,2])
        qp_pairs[j] = -1
        
    # Compute magnitudes of p and find the shortest r1
    for k in range(pnum):
        pmag[k] = sqrt(p[k,0] * p[k,0] + p[k,1] * p[k,1] + p[k,2] * p[k,2])
        if pmag[k] < r1:
            r1 = pmag[k]
    
    # Loop over all q
    for j in range(qnum):
        cos_theta_min = cos_theta_max
    
        # Loop over all p
        for k in range(pnum):
            
            # Compute cos theta between p and q
            cos_theta = (q[j,0] * p[k,0] + q[j,1] * p[k,1] + q[j,2] * p[k,2]) / (qmag[j] * pmag[k])
            
            # Select q-p pair with the smallest theta angle (largest cos theta)
            if cos_theta > cos_theta_min:
                cos_theta_min = cos_theta
                qp_pairs[j] = k
                
        # Check if the best p is already assigned to another q
        if qp_pairs[j] >= 0:
            for k in range(j):
                if qp_pairs[j] == qp_pairs[k]:
                    jrad = fabs(r1 - qmag[j])
                    krad = fabs(r1 - qmag[k])

                    # Remove the p-q pair that is farther from r1
                    if jrad < krad:
                        qp_pairs[k]=-1
                    else:
                        qp_pairs[j]=-1 
        
    # Construct reduced P, Q matrices from p-q pairs
    n = 0
    for j in range(qnum):
        if qp_pairs[j] >= 0:
            for x in range(3):
                Q[n,x] = q[j,x]
                P[n,x] = p[qp_pairs[j],x]
            n+=1
            
    # Return number of matching P-Q pairs 
    return n

@cython.boundscheck(False)
@cython.wraparound(False) 
cdef strain_c(const double[:,:,::1] G):
    """
    Computes the per-atom strain tensor from the lattice correspondence tensor G
    """
    # Parameter initialization 
    natoms = G.shape[0]
    ident = np.identity(3)
    cdef double [:,::1] identity = ident
    
    strain = np.empty((natoms, 3, 3))
    cdef double [:,:,::1] strain_view = strain
    cdef Py_ssize_t i, j, k
    
    # ε_jk = ((I_jk - G_jk) + (I_kj - G_kj)) / 2
    for i in range(natoms):
        for j in range(3):
            for k in range(3):
                strain_view[i,j,k] = ((identity[j,k] - G[i,j,k]) + (identity[k,j] - G[i,k,j])) / 2.
    return strain

@cython.boundscheck(False)
@cython.wraparound(False) 
cdef invariant1_c(const double [:,:,::1] strain):
    """
    Computes the per-atom first strain invariant from the strain tensor
    """
    # Parameter initialization 
    natoms = strain.shape[0]
    
    invar = np.empty(natoms)
    cdef double [::1] invar_view = invar
    cdef Py_ssize_t i
    
    # I1 = trace(ε_ij)
    for i in range(natoms):
        invar_view[i] = strain[i,0,0] + strain[i,1,1] + strain[i,2,2]
    
    return invar

@cython.boundscheck(False)
@cython.wraparound(False) 
cdef invariant2_c(const double [:,:,::1] strain):
    """
    Computes the per-atom second strain invariant from the strain tensor
    """
    # Parameter initialization 
    natoms = strain.shape[0]
    
    invar = np.empty(natoms)
    cdef double [::1] invar_view = invar
    cdef Py_ssize_t i
    
    # I2 = trace(ε_ij)
    for i in range(natoms):
        invar_view[i] = ( strain[i,0,0] * strain[i,1,1] 
                        + strain[i,0,0] * strain[i,2,2] 
                        + strain[i,1,1] * strain[i,2,2] 
                        - strain[i,0,1] * strain[i,1,0] 
                        - strain[i,0,2] * strain[i,2,0]
                        - strain[i,1,2] * strain[i,2,1])
    
    return invar

@cython.boundscheck(False)
@cython.wraparound(False) 
cdef invariant3_c(const double [:,:,::1] strain):
    """
    Computes the per-atom third strain invariant from the strain tensor
    """
    # Parameter initialization 
    natoms = strain.shape[0]
    
    invar = np.empty(natoms)
    cdef double [::1] invar_view = invar
    cdef Py_ssize_t i
    
    # I3 = det(ε_ij)
    for i in range(natoms):
        invar_view[i] = ( strain[i,0,0] * ( strain[i,1,1] * strain[i,2,2]
                                          - strain[i,1,2] * strain[i,2,1]) 
                        - strain[i,0,1] * ( strain[i,1,0] * strain[i,2,2]
                                          - strain[i,1,2] * strain[i,2,0]) 
                        + strain[i,0,2] * ( strain[i,1,0] * strain[i,2,1]
                                          - strain[i,1,1] * strain[i,2,0]))
    
    return invar

@cython.boundscheck(False)
@cython.wraparound(False) 
cdef rotation_c(const double[:,:,::1] G):
    """
    Computes the per-atom rotation tensor from the lattice correspondence tensor G
    """
    # Parameter initialization 
    natoms = G.shape[0]
    ident = np.identity(3)
    cdef double [:,::1] identity = ident
    
    rotation = np.empty((natoms, 3, 3))
    cdef double [:,:,::1] rot_view = rotation
    cdef Py_ssize_t i, j, k
    
    # rot_jk = ((I_jk - G_jk) - (I_kj - G_kj)) / 2
    for i in range(natoms):
        for j in range(3):
            for k in range(3):
                rot_view[i,j,k] = ((identity[j,k] - G[i,j,k]) - (identity[k,j] - G[i,k,j])) / 2.
    return rotation

@cython.boundscheck(False)
@cython.wraparound(False) 
cdef angularvelocity_c(const double [:,:,::1] rotation):
    """
    Computes the per-atom angular velocity from the rotation tensor
    """
    # Parameter initialization 
    natoms = rotation.shape[0]
    
    angvel = np.empty(natoms)
    cdef double [::1] angvel_view = angvel
    cdef Py_ssize_t i
    
    # a = sqrt(rot_12^2 + rot_13^2 + rot_23^2)
    for i in range(natoms):
        angvel_view[i] = sqrt( rotation[i,0,1] * rotation[i,0,1]
                             + rotation[i,0,2] * rotation[i,0,2]
                             + rotation[i,1,2] * rotation[i,1,2] )
    return angvel   



@cython.boundscheck(False)
@cython.wraparound(False) 
cdef dG_c(const double[:,:,::1] G, const long long[:,::1] nlist, double[:,:,::1] dG, Py_ssize_t i):
    """
    Computes the change in G between an atom and each of its neighbors
    """
    coord = nlist[i, 0]    
    cdef Py_ssize_t j, x, y
    
    # Calculate change in G between all neighbors
    for j in range(coord):
        for x in range(3):
            for y in range(3):
                dG[j, x, y] = G[nlist[i, j+1], x, y] - G[i, x, y]
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef nye_c(double [:,:,::1] gradG, double [:,:,::1] nye, Py_ssize_t i):
    """
    Computes the per-atom Nye tensor from the calculated grad(G).
    """
    nye[i,0,0] = gradG[1,0,2] - gradG[2,0,1]
    nye[i,0,1] = gradG[1,1,2] - gradG[2,1,1]
    nye[i,0,2] = gradG[1,2,2] - gradG[2,2,1]
    nye[i,1,0] = gradG[2,0,0] - gradG[0,0,2]
    nye[i,1,1] = gradG[2,1,0] - gradG[0,1,2]
    nye[i,1,2] = gradG[2,2,0] - gradG[0,2,2]
    nye[i,2,0] = gradG[0,0,1] - gradG[1,0,0]
    nye[i,2,1] = gradG[0,1,1] - gradG[1,1,0]
    nye[i,2,2] = gradG[0,2,1] - gradG[1,2,0]