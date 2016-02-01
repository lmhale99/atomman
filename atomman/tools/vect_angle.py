import numpy as np

def vect_angle(vect1, vect2):
    """Returns the angle between two vectors in degrees."""
    vect1 = np.asarray(vect1)
    vect2 = np.asarray(vect2)
    u_vect1 = vect1 / np.linalg.norm(vect1)
    u_vect2 = vect2 / np.linalg.norm(vect2)
    cosine = np.dot(u_vect1, u_vect2)
    try:
        cosine[cosine < -1] = -1
        cosine[cosine > 1] = 1
    except:
        if cosine < -1: cosine = -1
        elif cosine > 1: cosine = 1
    
    return 180 * np.arccos(cosine) / np.pi