# coding: utf-8

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

__all__ = ['plane3to4', 'plane4to3', 'vector3to4', 'vector4to3',
           'plane_crystal_to_cartesian',
           'vector_crystal_to_cartesian',
           'vector_primitive_to_conventional', 
           'vector_conventional_to_primitive']

def plane3to4(indices: npt.ArrayLike) -> np.ndarray:
    """
    Converts 3-term Miller (hkl) plane indices to 4-term hexagonal (hkil)
    Miller-Bravias indices.
    
    Parameters
    ----------
    indices : array-like object
        (..., 3) array of Miller crystallographic indices.
   
    Returns
    -------
    numpy.ndarray
        (..., 4) array of Miller-Bravais crystallographic indices.
    
    Raises
    ------
    ValueError
        If indices dimensions are not (..., 3).
    """
    
    # Verify compatible indices
    indices = np.asarray(indices)
    if indices.shape[-1] != 3:
        raise ValueError('Invalid index dimensions')
    
    newindices = np.empty(indices.shape[:-1] + (4,))
    newindices[..., 0] = indices[..., 0]
    newindices[..., 1] = indices[..., 1]
    newindices[..., 2] = -(newindices[..., 0] + newindices[..., 1])
    newindices[..., 3] = indices[..., 2]
    
    return newindices

def plane4to3(indices: npt.ArrayLike) -> np.ndarray:
    """
    Converts 4-term hexagonal Miller-Bravias (hkil) plane indices to 3-term 
    Miller (hkl) indices.
    
    Parameters
    ----------
    indices : array-like object
        (..., 4) array of Miller-Bravais crystallographic indices.
   
    Returns
    -------
    numpy.ndarray
        (..., 3) array of Miller crystallographic indices.
   
    Raises
    ------    
    ValueError
        If indices dimensions are not (..., 4), or if h+k+i != 0.
    """
    
    # Verify compatible indices
    indices = np.asarray(indices)
    if indices.shape[-1] != 4:
        raise ValueError('Invalid index dimensions')
    if not np.allclose(indices[...,:3].sum(axis=-1), 0.0):
        raise ValueError('Invalid indices: h+k+i != 0')
    
    newindices = np.empty(indices.shape[:-1] + (3,))
    newindices[..., 0] = indices[..., 0]
    newindices[..., 1] = indices[..., 1]
    newindices[..., 2] = indices[..., 3]
    
    return newindices

def vector3to4(indices: npt.ArrayLike) -> np.ndarray:
    """
    Converts 3-term Miller [uvw] vector indices to 4-term hexagonal [uvtw]
    Miller-Bravias indices. 
    
    Parameters
    ----------
    indices : array-like object
        (..., 3) array of Miller crystallographic indices.
   
    Returns
    -------
    numpy.ndarray
        (..., 4) array of Miller-Bravais crystallographic indices.
    
    Raises
    ------
    ValueError
        If indices dimensions are not (..., 3).
    """
    # Verify compatible indices
    indices = np.asarray(indices)
    if indices.shape[-1] != 3:
        raise ValueError('Invalid index dimensions')
    
    # Transform 
    newindices = np.empty(indices.shape[:-1] + (4,))
    newindices[..., 0] = (2 * indices[..., 0] - indices[..., 1]) / 3
    newindices[..., 1] = (2 * indices[..., 1] - indices[..., 0]) / 3
    newindices[..., 2] = -(newindices[..., 0] + newindices[..., 1])
    newindices[..., 3] = indices[..., 2] 
    
    return newindices
    
def vector4to3(indices: npt.ArrayLike) -> np.ndarray:
    """
    Converts 4-term hexagonal Miller-Bravias [uvtw] vector indices to 3-term 
    Miller [uvw] indices.
    
    Parameters
    ----------
    indices : array-like object
        (..., 4) array of Miller-Bravais crystallographic indices.
   
    Returns
    -------
    numpy.ndarray
        (..., 3) array of Miller crystallographic indices.
   
    Raises
    ------    
    ValueError
        If indices dimensions are not (..., 4), or if u+v+t != 0.
    """
    # Verify compatible indices
    indices = np.asarray(indices)
    if indices.shape[-1] != 4:
        raise ValueError('Invalid index dimensions')
    if not np.allclose(indices[...,:3].sum(axis=-1), 0.0):
        raise ValueError('Invalid indices: u+v+t != 0')
    
    # Transform
    newindices = np.empty(indices.shape[:-1] + (3,))
    newindices[..., 0] = 2 * indices[..., 0] + indices[..., 1]
    newindices[..., 1] = 2 * indices[..., 1] + indices[..., 0]
    newindices[..., 2] = indices[..., 3] 
    
    return newindices

def vector_crystal_to_cartesian(indices: npt.ArrayLike,
                                box) -> np.ndarray:
    """
    Converts crystal indices to Cartesian vectors relative
    to a given lattice box. 
    
    Parameters
    ----------
    indices : array-like object
        (..., 3) array of [uvw] Miller crystallographic indices or 
        (..., 4) array of [uvtw] Miller-Bravais crystallographic indices.
    box : atomman.Box
        Box that defines the lattice cell vectors to use. 
   
    Returns
    -------
    numpy.ndarray
        (..., 3) array of Cartesian vectors.
        
    Raises
    ------
    ValueError
        If indices dimensions are not (..., 3) or (..., 4), or if
        hexagonal indices given with non-hexagonal box.
    """
    indices = np.array(indices)
    
    # Convert 4-term Miller-Bravais to standard 3-term indices
    if indices.shape[-1] == 4: 
        if not box.ishexagonal():
            raise ValueError('Hexagonal indices given with non-hexagonal box')
        indices = vector4to3(indices)
    
    if indices.shape[-1] != 3:
        raise ValueError('Invalid index dimensions')

    return indices.dot(box.vects)
    
def vector_primitive_to_conventional(indices: npt.ArrayLike,
                                     setting: str = 'p') -> np.ndarray:
    """
    Converts crystal indices relative to a primitive cell 
    to indices relative to a conventional cell in a specified setting.
    
    Parameters
    ----------
    indices : array-like object
        (..., 3) array of [uvw] Miller crystallographic indices
        relative to the primitive cell
    setting : str
        Specifies the conventional cell setting: 'p' for primitive,
        'a', 'b', 'c' for side-centered, 'i' for body-centered, and
        'f' for face-centered.
   
    Returns
    -------
    numpy.ndarray
        (..., 3) array of [uvw] Miller crystallographic indices
        relative to the conventional cell
        
    Raises
    ------
    ValueError
        If indices dimensions are not (..., 3) or if an unknown setting
        value is given.
    """

    # Verify compatible indices
    indices = np.asarray(indices)
    if indices.shape[-1] != 3:
        raise ValueError('Invalid index dimensions')

    lattice_vectors = {}
    lattice_vectors['p'] = np.array([[ 1.0, 0.0, 0.0],
                                     [ 0.0, 1.0, 0.0],
                                     [ 0.0, 0.0, 1.0]])
    lattice_vectors['a'] = np.array([[ 1.0, 0.0, 0.0],
                                     [ 0.0, 0.5, 0.5],
                                     [ 0.0,-0.5, 0.5]])
    lattice_vectors['b'] = np.array([[ 0.5, 0.0, 0.5],
                                     [ 0.0, 1.0, 0.0],
                                     [-0.5, 0.0, 0.5]])
    lattice_vectors['c'] = np.array([[ 0.5, 0.5, 0.0],
                                     [-0.5, 0.5, 0.0],
                                     [ 0.0, 0.0, 1.0]])
    lattice_vectors['i'] = np.array([[ 0.5, 0.5, 0.5],
                                     [-0.5, 0.5,-0.5],
                                     [-0.5,-0.5, 0.5]])
    lattice_vectors['f'] = np.array([[ 0.5, 0.5, 0.0],
                                     [ 0.0, 0.5, 0.5],
                                     [ 0.5, 0.0, 0.5]])
    
    try:
        lat = lattice_vectors[setting]
    except:
        raise ValueError('Unknown lattice setting. Allowed values are: p, a, b, c, i, and f')
    
    return indices.dot(lat)

def vector_conventional_to_primitive(indices: npt.ArrayLike,
                                     setting: str = 'p') -> np.ndarray:
    """
    Converts crystal indices relative to a conventional cell 
    in a specified setting to indices relative to a primitive cell.
    
    Parameters
    ----------
    indices : array-like object
        (..., 3) array of [uvw] Miller crystallographic indices
        relative to the conventional cell
    setting : str
        Specifies the conventional cell setting: 'p' for primitive,
        'a', 'b', 'c' for side-centered, 'i' for body-centered, and
        'f' for face-centered.
   
    Returns
    -------
    numpy.ndarray
        (..., 3) array of [uvw] Miller crystallographic indices
        relative to the primitive cell
        
    Raises
    ------
    ValueError
        If indices dimensions are not (..., 3) or if an unknown setting
        value is given.
    """
    # Verify compatible indices
    indices = np.asarray(indices)
    if indices.shape[-1] != 3:
        raise ValueError('Invalid index dimensions')
    
    lattice_vectors = {}
    lattice_vectors['p'] = np.array([[ 1.0, 0.0, 0.0],
                                     [ 0.0, 1.0, 0.0],
                                     [ 0.0, 0.0, 1.0]])
    lattice_vectors['a'] = np.array([[ 1.0, 0.0, 0.0],
                                     [ 0.0, 1.0,-1.0],
                                     [ 0.0, 1.0, 1.0]])
    lattice_vectors['b'] = np.array([[ 1.0, 0.0,-1.0],
                                     [ 0.0, 1.0, 0.0],
                                     [ 1.0, 0.0, 1.0]])
    lattice_vectors['c'] = np.array([[ 1.0,-1.0, 0.0],
                                     [ 1.0, 1.0, 0.0],
                                     [ 0.0, 0.0, 1.0]])
    lattice_vectors['i'] = np.array([[ 0.0,-1.0,-1.0],
                                     [ 1.0, 1.0, 0.0],
                                     [ 1.0, 0.0, 1.0]])
    lattice_vectors['f'] = np.array([[ 1.0,-1.0, 1.0],
                                     [ 1.0, 1.0,-1.0],
                                     [-1.0, 1.0, 1.0]])
    
    try:
        lat = lattice_vectors[setting]
    except:
        raise ValueError('Unknown lattice setting. Allowed values are: p, a, b, c, i, and f')
    
    return indices.dot(lat)

def plane_crystal_to_cartesian(indices: npt.ArrayLike,
                               box) -> np.ndarray:
    """
    Converts crystal planar indices to Cartesian plane normal vectors relative
    to a given lattice box.  Note: the algorithm used requires that the planar
    indices be integers.
    
    Parameters
    ----------
    indices : array-like object
        (..., 3) array of [hkl] Miller crystallographic indices or 
        (..., 4) array of [hkil] Miller-Bravais crystallographic indices.
    box : atomman.Box
        Box that defines the lattice cell vectors to use. 
   
    Returns
    -------
    numpy.ndarray
        (..., 3) array of Cartesian vectors corresponding to plane normals.
        
    Raises
    ------
    ValueError
        If indices dimensions are not (..., 3) or (..., 4), or if
        hexagonal indices given with non-hexagonal box.
    """
    # Check indices values
    indices = np.asarray(indices)
    
    # Convert hkil to hkl
    if indices.shape[-1] == 4:
        if box.ishexagonal():
            indices = plane4to3(indices)
        else:
            raise ValueError('Hexagonal indices given with non-hexagonal box')
    if indices.shape[-1] != 3:
        raise ValueError('Invalid index dimensions')
    
    # Verify values are int
    if np.allclose(indices, np.asarray(indices, dtype=int)):
        indices = np.asarray(indices, dtype=int)
    else:
        raise ValueError('Indices must be integers')
    
    
    def plane_cryst_2_cart(indices, box):
        """
        Finds the Cartesian plane normal for a single 3 index Miller crystal
        plane and a box.        

        Parameters
        ----------
        indices : numpy.ndarray
            The 3 index Miller crystal plane
        box : atomman.Box
            The unit cell box.

        Returns
        -------
        numpy.ndarray
            The Cartesian plane normal vector.
        """
        # Find two in-plane box vectors
        if indices[0] != 0:
            if indices[1] != 0:
                if indices[2] != 0:
                    # indices solution
                    m = np.lcm.reduce([indices[0], indices[1], indices[2]])
                    s = np.sign(indices[0] * indices[1] * indices[2])
                    a_uvw = np.array([-m / indices[0], m / indices[1], 0], dtype=int)
                    b_uvw = np.array([-m / indices[0], 0, m / indices[2]], dtype=int)
                else:
                    # hk0 solution
                    m = np.lcm(indices[0], indices[1])
                    s = np.sign(indices[0] * indices[1])
                    a_uvw = np.array([-m / indices[0], m / indices[1], 0], dtype=int)
                    b_uvw = np.array([0, 0, 1], dtype=int)
            else:
                if indices[2] != 0:
                    # h0l solution
                    m = np.lcm(indices[0], indices[2])
                    s = np.sign(indices[0] * indices[2])
                    a_uvw = np.array([m / indices[0], 0, -m / indices[2]], dtype=int)
                    b_uvw = np.array([0, 1, 0], dtype=int)
                else:
                    # h00 solution
                    m = 1
                    s = np.sign(indices[0])
                    a_uvw = np.array([0, 1, 0], dtype=int)
                    b_uvw = np.array([0, 0, 1], dtype=int)
        elif indices[1] != 0:
            if indices[2] != 0:
                # 0kl solution
                m = np.lcm(indices[1], indices[2]) 
                s = np.sign(indices[1] * indices[2])
                a_uvw = np.array([0, -m / indices[1], m / indices[2]], dtype=int)
                b_uvw = np.array([1, 0, 0], dtype=int)
            else:    
                # 0k0 solution
                m = 1
                s = np.sign(indices[1])
                a_uvw = np.array([0, 0, 1], dtype=int)
                b_uvw = np.array([1, 0, 0], dtype=int)
        elif indices[2] != 0:
            # 00l solution
            m = 1
            s = np.sign(indices[2])
            a_uvw = np.array([1, 0, 0], dtype=int)
            b_uvw = np.array([0, 1, 0], dtype=int)
        else:
            raise ValueError('indices cannot be all zeros')
        
        # Compute Cartesian plane normal
        planenormal = s * np.cross(a_uvw.dot(box.vects), b_uvw.dot(box.vects))

        # Return the unit vector normal
        return planenormal / np.linalg.norm(planenormal)

    # Apply plane_cryst_2_cart to each given set of indices
    return np.apply_along_axis(plane_cryst_2_cart, -1, indices, box)