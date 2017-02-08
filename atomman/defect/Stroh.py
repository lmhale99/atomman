import numpy as np
from atomman.tools import ElasticConstants, axes_check
from copy import deepcopy

class Stroh(object):
    """
    Class for solving the Eshelby anisotropic solution for a perfectly stright
    dislocation or crack using the Stroh method.
    
    Class Methods:
    solve(C, b, axes=None, tol=1e-8) -- solves the elasticity problem for given 
                                        elastic constants, Burgers vector, 
                                        (axes, and rounding tolerance)
    displacement(x, y) -- compute the position-dependent 
                                anisotropic displacement
    stress(x, y) -- compute the position-dependent anisotropic stress 
                          state
    
    Class Properties (These values cannot be set directly):
    preln -- the pre-ln elastic energy factor. 
    p -- the eigenvalues of the solution.
    A -- the A eigenvectors of the solution.
    L -- the L eigenvectors of the solution.
    k -- normalization term.
    """
    
    def __init__(self, C, burgers, axes=None, tol=1e-8): 
        """
        Initialize an instance of the Stroh class and calls solve().
        
        Required Arguments:
        C -- an instance of ElasticConstants.
        b -- a numpy array representing the Burgers vector.
        
        Keyword Arguments:
        axes -- rotational axes of the system. If given, then C and b will be 
                transformed.
        tol -- tolerance parameter used to round off near-zero values.
        """

        self.solve(C, burgers, axes=axes, tol=tol)

    def solve(self, C, burgers, axes=None, tol=1e-8):
        """Performs the Stroh method to obtain solution terms p, A, L, and k.
        
        Required Arguments:
        C -- an instance of ElasticConstants.
        b -- a numpy array representing the Burgers vector.
        
        Keyword Arguments:
        axes -- rotational axes of the system. If given, then C and b will be 
                transformed.
        tol -- tolerance parameter used to round off near-zero values.
        """
        burgers = np.asarray(burgers, dtype='float64')
        if axes is not None:
            T = axes_check(axes)
            burgers = T.dot(burgers)
            C = C.transform(axes)
        Cijkl = C.Cijkl   
        
        m = np.array([1.0, 0.0, 0.0])
        n = np.array([0.0, 1.0, 0.0])
        
        #Matrixes of Cijkl constants used to construct N
        mm = np.einsum('i,ijkl,l', m, Cijkl, m)
        mn = np.einsum('i,ijkl,l', m, Cijkl, n)
        nm = np.einsum('i,ijkl,l', n, Cijkl, m)
        nn = np.einsum('i,ijkl,l', n, Cijkl, n)
        
        #The four 3x3 matrixes that represent the quadrants of N
        NB = -np.linalg.inv(nn)
        NA = NB.dot(nm)
        NC = mn.dot(NA) + mm
        ND = mn.dot(NB)
        
        #N is the 6x6 array, where the eigenvalues are the roots p
        #and the eigenvectors give A and L
        N =  np.array(np.vstack((np.hstack((NA, NB)), np.hstack((NC, ND)))))
        
        #Calculate the eigenvectors and eigenvalues
        eig = np.linalg.eig(N)
        p = eig[0]
        eigvec = np.transpose(eig[1])

        #separate the eigenvectors into A and L
        A = np.array([eigvec[0,:3], eigvec[1,:3], eigvec[2,:3], eigvec[3,:3], eigvec[4,:3], eigvec[5,:3]])
        L = np.array([eigvec[0,3:], eigvec[1,3:], eigvec[2,3:], eigvec[3,3:], eigvec[4,3:], eigvec[5,3:]])
        
        #calculate k
        k = 1. / (2. * np.einsum('si,si->s', A, L))
        
        #Calculation verification checks
        #try:
        #    assert np.allclose(np.einsum('s,si,sj->ij', k,A,L), np.identity(3, dtype='complex128'), atol=tol)
        #    assert np.allclose(np.einsum('s,si,sj->ij', k,A,A), np.zeros((3,3), dtype='complex128'), atol=tol)
        #    assert np.allclose(np.einsum('s,si,sj->ij', k,L,L), np.zeros((3,3), dtype='complex128'), atol=tol)
        #    assert np.allclose(np.einsum('s,t,si,ti->st', k**.5,k**.5,A,L) + np.einsum('s,t,ti,si->st', k**.5,k**.5,A,L),
        #                                 np.identity(6, dtype='complex128'), atol = tol)
        #except:
        #    raise ValueError('Stroh checks failed!')
    
        #assign property values
        self.__burgers = burgers
        self.__Cijkl = Cijkl
        self.__tol = tol
        self.__p = p
        self.__A = A
        self.__L = L
        self.__k = k
        self.__m = m
        self.__n = n
        
    @property
    def p(self):
        """p eigenvalues"""
        return deepcopy(self.__p)
        
    @property
    def A(self):
        """A eigenvectors"""
        return deepcopy(self.__A)
        
    @property
    def L(self):
        """L eigenvectors"""
        return deepcopy(self.__L)
        
    @property
    def k(self):
        """k normalization factors"""
        return deepcopy(self.__k)

    @property
    def K_coeff(self):
        """The energy coefficient"""
        
        return self.__burgers.dot(self.K_tensor.dot(self.__burgers)) / self.__burgers.dot(self.__burgers)        
        
    @property
    def K_tensor(self):
        """The energy coefficient tensor"""
        
        ii = np.array([1.j])
        updn = np.array([1, -1, 1, -1, 1, -1])
        
        K = ii * np.einsum('s,s,si,sj->ij', updn, self.k, self.L, self.L)
        K = np.real_if_close(K, tol=self.__tol)
        
        K[np.isclose(K / K.max(), 0.0, atol=self.__tol)] = 0.0
        
        return K
        
    @property
    def preln(self):
        """The pre-ln strain energy factor"""
        
        ii = np.array([1.j])
        updn = np.array([1, -1, 1, -1, 1, -1])
        
        factor = ii * np.einsum('s,s,si,sj->ij', updn, self.k, self.L, self.L) / (4 * np.pi)
        factor = np.real_if_close(factor, tol=self.__tol)
        
        return self.__burgers.dot(factor.dot(self.__burgers))

    def displacement(self, pos):
        """
        Compute the position-dependent anisotropic displacement. 
        
        The argument pos can either be a single 3D position,
        or a list/array of 3D positions.
        """
        
        ii = np.array([1.j])
        updn = np.array([1, -1, 1, -1, 1, -1])
        eta = self.eta(pos)
        
        disp = 1 / (2 * np.pi * ii) * np.einsum('a,a,ai,a,na->ni', updn, self.k, self.A, self.L.dot(self.__burgers), np.log(eta))
        disp = np.real_if_close(disp, tol=self.__tol)
        
        if np.asarray(pos).ndim == 1:
            return disp[0]
        else:
            return disp

    def stress(self, pos):
        """
        Compute the position-dependent anisotropic stress state. 
        
        The argument pos can either be a single 3D position,
        or a list/array of 3D positions.
        """
        
        ii = np.array([1.j])
        updn = np.array([1, -1, 1, -1, 1, -1])
        eta = self.eta(pos)
        
        mpn = self.__m + np.outer(self.p, self.__n)
        stress = 1 / (2 * np.pi * ii) * np.einsum('a,a,ijkl,al,ak,a,na->nij', updn, self.k, self.__Cijkl, mpn, self.A, self.L.dot(self.__burgers), 1/eta)
        stress = np.real_if_close(stress, tol=self.__tol)
        
        if np.asarray(pos).ndim == 1:
            return stress[0]
        else:
            return stress

    def eta(self, pos):
        """Compute the position dependent eta factor"""
        
        x = np.dot(pos, self.__m)
        y = np.dot(pos, self.__n)
        
        return (x + np.outer(self.p, y)).T