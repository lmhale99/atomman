# Standard Python libraries
from __future__ import (absolute_import, print_function,
                        division, unicode_literals)

# http://www.numpy.org/
import numpy as np

from .crystalsystem import ishexagonal

__all__ = ['vector3to4', 'vector4to3', 'vector_crystal_to_cartesian',
           'vector_primitive_to_conventional', 
           'vector_conventional_to_primitive']

def vector3to4(indices):
    """
    Converts 3-term Miller [uvw] indices to 4-term hexagonal [uvtw]
    Miller-Bravias indices. 
    
    Parameters
    ----------
    indices : np.ndarray
        (..., 3) array of Miller crystallographic indices.
   
    Returns
    -------
    np.ndarray
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
    
def vector4to3(indices):
    """
    Converts 4-term hexagonal Miller-Bravias [uvtw] indices to 3-term 
    Miller [uvw] indices. Note that vectors will be normalized to 
    smallest integer representations.
    
    Parameters
    ----------
    indices : np.ndarray of int
        (..., 4) array of Miller-Bravais crystallographic indices.
   
    Returns
    -------
    np.ndarray of int
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
    
    # Transform
    newindices = np.empty(indices.shape[:-1] + (3,))
    newindices[..., 0] = 2 * indices[..., 0] + indices[..., 1]
    newindices[..., 1] = 2 * indices[..., 1] + indices[..., 0]
    newindices[..., 2] = indices[..., 3] 
    
    return newindices

def vector_crystal_to_cartesian(indices, box):
    """
    Converts crystal indices to Cartesian vectors relative
    to a given lattice box. 
    
    Parameters
    ----------
    indices : np.ndarray
        (..., 3) array of [uvw] Miller crystallographic indices or 
        (..., 4) array of [uvtw] Miller-Bravais crystallographic indices.
    box : atomman.Box
        Box that defines the lattice cell vectors to use. 
   
    Returns
    -------
    np.ndarray of float
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
        if not ishexagonal(box):
            raise ValueError('Hexagonal indices given with non-hexagonal box')
        indices = vector4to3(indices)
    
    if indices.shape[-1] != 3:
        raise ValueError('Invalid index dimensions')

    return (np.outer(indices[...,0], box.vects[0]) 
          + np.outer(indices[...,1], box.vects[1]) 
          + np.outer(indices[...,2], box.vects[2])).reshape(indices.shape)
    
def vector_primitive_to_conventional(indices, setting='p'):
    """
    Converts crystal indices relative to a primitive cell 
    to indices relative to a conventional cell in a specified setting.
    
    Parameters
    ----------
    indices : np.ndarray
        (..., 3) array of [uvw] Miller crystallographic indices
        relative to the primitive cell
    setting : str
        Specifies the conventional cell setting: 'p' for primitive,
        'a', 'b', 'c' for side-centered, 'i' for body-centered, and
        'f' for face-centered.
   
    Returns
    -------
    np.ndarray of float
        (..., 3) array of [uvw] Miller crystallographic indices
        relative to the conventional cell
        
    Raises
    ------
    ValuenError
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
                                     [ 0.5, 0.0, 0.5],
                                     [ 0.0, 0.5, 0.5]])
    
    try:
        lat = lattice_vectors[setting]
    except:
        raise ValueError('Unknown lattice setting. Allowed values are: p, a, b, c, i, and f')
    
    return np.inner(indices, lat.T)


def vector_conventional_to_primitive(indices, setting='p'):
    """
    Converts crystal indices relative to a conventional cell 
    in a specified setting to indices relative to a primitive cell.
    
    Parameters
    ----------
    indices : np.ndarray
        (..., 3) array of [uvw] Miller crystallographic indices
        relative to the conventional cell
    setting : str
        Specifies the conventional cell setting: 'p' for primitive,
        'a', 'b', 'c' for side-centered, 'i' for body-centered, and
        'f' for face-centered.
   
    Returns
    -------
    np.ndarray of float
        (..., 3) array of [uvw] Miller crystallographic indices
        relative to the primitive cell
        
    Raises
    ------
    ValuenError
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
    lattice_vectors['f'] = np.array([[ 1.0, 1.0,-1.0],
                                     [ 1.0,-1.0, 1.0],
                                     [-1.0, 1.0, 1.0]])
    
    try:
        lat = lattice_vectors[setting]
    except:
        raise ValueError('Unknown lattice setting. Allowed values are: p, a, b, c, i, and f')
    
    return np.inner(indices, lat.T)