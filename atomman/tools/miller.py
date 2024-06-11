# coding: utf-8

# http://www.numpy.org/
import numpy as np
import numpy.typing as npt

__all__ = ['plane3to4', 'plane4to3', 'vector3to4', 'vector4to3',
           'plane_crystal_to_cartesian',
           'vector_crystal_to_cartesian',
           'vector_primitive_to_conventional', 
           'vector_conventional_to_primitive',
           'fromstring', 'tostring',
           'all_indices', 'reduce_indices']

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
        'a', 'b', 'c' for side-centered, 'i' for body-centered,
        'f' for face-centered, and 't' for trigonal systems.

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
    lattice_vectors['p'] = np.array([[  1.0,  0.0,  0.0],
                                     [  0.0,  1.0,  0.0],
                                     [  0.0,  0.0,  1.0]])
    
    lattice_vectors['a'] = np.array([[  1.0,  0.0,  0.0],
                                     [  0.0,  0.5,  0.5],
                                     [  0.0, -0.5,  0.5]])
    
    lattice_vectors['b'] = np.array([[  0.5,  0.0,  0.5],
                                     [  0.0,  1.0,  0.0],
                                     [ -0.5,  0.0,  0.5]])
    
    lattice_vectors['c'] = np.array([[  0.5,  0.5,  0.0],
                                     [ -0.5,  0.5,  0.0],
                                     [  0.0,  0.0,  1.0]])
    
    lattice_vectors['i'] = np.array([[  0.5,  0.5,  0.5],
                                     [ -0.5,  0.5, -0.5],
                                     [ -0.5, -0.5,  0.5]])
    
    lattice_vectors['f'] = np.array([[  0.5,  0.5,  0.0],
                                     [  0.0,  0.5,  0.5],
                                     [  0.5,  0.0,  0.5]])
    
    lattice_vectors['t1'] = np.array([[  2.0,  1.0,  1.0],
                                      [ -1.0,  1.0,  1.0],
                                      [ -1.0, -2.0,  1.0]]) / 3.
    
    lattice_vectors['t2'] = np.array([[ -2.0, -1.0,  1.0],
                                      [  1.0, -1.0,  1.0],
                                      [  1.0,  2.0,  1.0]]) / 3.

    try:
        lat = lattice_vectors[setting]
    except:
        raise ValueError('Unknown lattice setting. Allowed values are: p, a, b, c, i, t and f')

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
        'a', 'b', 'c' for side-centered, 'i' for body-centered,
        'f' for face-centered, and 't1' or 't2' for trigonal systems.

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
    lattice_vectors['p'] = np.array([[  1.0,  0.0,  0.0],
                                     [  0.0,  1.0,  0.0],
                                     [  0.0,  0.0,  1.0]])
    
    lattice_vectors['a'] = np.array([[  1.0,  0.0,  0.0],
                                     [  0.0,  1.0, -1.0],
                                     [  0.0,  1.0,  1.0]])
    
    lattice_vectors['b'] = np.array([[  1.0,  0.0, -1.0],
                                     [  0.0,  1.0,  0.0],
                                     [  1.0,  0.0,  1.0]])
    
    lattice_vectors['c'] = np.array([[  1.0, -1.0,  0.0],
                                     [  1.0,  1.0,  0.0],
                                     [  0.0,  0.0,  1.0]])
    
    lattice_vectors['i'] = np.array([[  0.0, -1.0, -1.0],
                                     [  1.0,  1.0,  0.0],
                                     [  1.0,  0.0,  1.0]])
    
    lattice_vectors['f'] = np.array([[  1.0, -1.0,  1.0],
                                     [  1.0,  1.0, -1.0],
                                     [ -1.0,  1.0,  1.0]])
    
    lattice_vectors['t1'] = np.array([[  1.0, -1.0,  0.0],
                                      [  0.0,  1.0, -1.0],
                                      [  1.0,  1.0,  1.0]])
    
    lattice_vectors['t2'] = np.array([[ -1.0,  1.0,  0.0],
                                      [  0.0, -1.0,  1.0],
                                      [  1.0,  1.0,  1.0]])

    try:
        lat = lattice_vectors[setting]
    except:
        raise ValueError('Unknown lattice setting. Allowed values are: p, a, b, c, i, t1, t2, and f')

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

def fromstring(value):
    """
    Reads Miller and Miller-Bravais vector/plane strings and returns a corresponding numpy
    array.  Values should be space-delimited and surrounded by angle brackets. The
    values in the brackets should be integers. Fractional arrays can be defined by placing
    a fraction before the array.
    
    Examples:
    [1 0 0]
    1/2 [1 1 0]
    [0 0 0 1]
    1/3 [1 1 -2 0]
    
    Parameters
    ----------
    value : str
        The Miller(-Bravais) vector/plane string to interpret.
        
    Returns
    -------
    numpy.NDArray
        The vector/plane as an array of floats.
    """
    # Find the str indices of any brackets
    if value.find('[') > -1:
        openindex = value.index('[')
        closeindex = value.index(']')
    elif value.find('(') > -1:
        openindex = value.index('(')
        closeindex = value.index(')')
    elif value.find('<') > -1:
        openindex = value.index('<')
        closeindex = value.index('>')
    elif value.find('{') > -1:
        openindex = value.index('{')
        closeindex = value.index('}')
    else:
        openindex = -1
        
    # Get leading fraction
    if openindex > 0:
        fraction = value[:openindex]
        terms = fraction.split('/')
        assert len(terms) == 2, 'fraction can only have one /'
        fraction = float(terms[0]) / float(terms[1])
    else:
        fraction = 1

     # Legacy reader for just numbers no brackets
    if openindex == -1:
        array = np.fromstring(value, dtype=float, sep=' ')

    # Convert the value in brackets into an array
    else:
        array = np.fromstring(value[openindex+1: closeindex], dtype=float, sep=' ')
    assert array.shape in [(3,), (4,)], 'array must have 3 or 4 indices'

    # Apply the fraction and return
    return fraction * array

def tostring(array, bracket='[]'):
    """
    Converts a 3-index Miller or 4-index Miller-Bravais vector/plane into a
    representative string.  This is the inverse of the fromstring operation.
    
    The values in the allowed vectors are limited to what is generally allowed
    for full lattice vectors/planes:
    - For planes, all values must be integers
    - For 3-index vectors, all values must be multiples of 1/2
    - For 4-index vectors, all values must be multiples of 1/3

    Note that this is a simple conversion of the numbers and does not actively
    check that the fractional vectors correspond to actual lattice vectors.
    Such checks would require using the other tools in the Miller module that
    convert between different representations of the vectors.

    Parameters
    ----------
    array : array-like object
        The 3-index Miller or 4-index Miller-Bravais vector/plane given as an
        array or list.
    bracket : str, optional
        The bracket style to use for surrounding the indices in the string.
        '[]' and '<>' indicate vectors, and '()' and '{}' indicate planes.

    Returns
    -------
    str
        The string representation of the vector/plane given as an optional
        fraction followed by bracketed space-delimited integers.
    """

    # Check dimensions
    array = np.asarray(array)
    if array.shape != (3,) and array.shape != (4,):
        raise ValueError('array must be 1D with 3 or 4 terms, i.e. a shape of (3,) or (4,)')

    # Handle brackets
    if bracket in ['[]', '<>']:
        vector = True
    elif bracket in ['()', '{}']:
        vector = False
    else:
        raise ValueError("Invalid bracket style: allowed styles are '[]', '<>', '()', and '{}'")

    # Convert float to int values
    if np.issubdtype(array.dtype, float):
        
        # Check if all values are ints (works for vectors and planes)
        intarray = np.asarray(np.around(array), dtype=int)
        if np.allclose(array, intarray):
            array = intarray
            frac = ''

        # Throw error for non-integer plane indices
        elif vector is False:
            raise ValueError('plane indices must be integers')

        # Check if 3-index vector values are multiples of 1/2
        elif array.shape == (3,):
            int2array = np.asarray(np.around(2 * array), dtype=int)
            if np.allclose(2 * array, int2array):
                array = int2array
                frac = '1/2 '
            else:
                raise ValueError('3-index vector values must be multiples of 1/2')

        # Check if 4-index vector values are multiples of 1/3
        else:
            int3array = np.asarray(np.around(3 * array), dtype=int)
            if np.allclose(3 * array, int3array):
                array = int3array
                frac = '1/3 '
            else:
                raise ValueError('4-index vector values must be multiples of 1/3')

    # Check that any non-float values are (converted to) ints
    else:
        array = np.asarray(array, int)
        frac = ''

    # Generate the string
    if array.shape == (3,):
        return f'{frac}{bracket[0]}{array[0]} {array[1]} {array[2]}{bracket[1]}'
    else:
        return f'{frac}{bracket[0]}{array[0]} {array[1]} {array[2]} {array[3]}{bracket[1]}'

def all_indices(maxindex: int = 10,
                reduce: bool = False) -> np.ndarray:
    """
    Generates an array containing all Miller indices where each individual
    index value is an integer independently ranging from -maxindex to maxindex. 

    Parameters
    ----------
    maxindex : int, optional
        The maximum absolute index to use: u,v,w (or h,k,l) can all
        independently vary from -maxindex to maxindex.  Default value is 10.
    reduce : bool, optional
        Setting this to True will return only the indices sets that are unique
        after reducing to their smallest integer values.  Default value of
        False will return all indices even those that are multiples of others.

    Returns 
    -------
    indices : numpy.NDArray
        All integer Miller [uvw] or (hkl) indices within the maxindex range.
    """
    # Individual indices independently range -maxindex to maxindex
    i = np.arange(-maxindex, maxindex+1)
    
    # Build grid, combine and transform
    u, v, w = np.meshgrid(i, i, i)
    indices = np.vstack([u.flat, v.flat, w.flat]).T
    
    # Remove [0,0,0]
    indices = indices[np.abs(indices).sum(axis=1) != 0]

    # Reduce and remove duplicates
    if reduce:
        indices = reduce_indices(indices)
        indices = np.unique(indices, axis=0)

    return indices

def reduce_indices(indices: npt.ArrayLike) -> np.array:
    """
    Given an array of one or more Miller(-Bravais) indices reduce them to the
    smallest integer representation.
    
    Parameters
    ----------
    indices : Array-like object
        An array of ints with shapes (...,3) or (...,4) that represent one
        or more Miller(-Bravais) crystal vectors.

    Returns
    -------
    numpy.ndarray
        The reduced indices
    """
    indices = np.asarray(indices)
    if indices.shape[-1] != 3 and indices.shape[-1] != 4:
        raise ValueError('Invalid indices dimensions')
    
    # Get gcd of each vector
    n = np.gcd.reduce(indices, axis=-1)

    # Divide each vector by the corresponding gcd
    red_indices = (indices.T // n).T

    return red_indices

