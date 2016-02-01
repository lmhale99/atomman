import numpy as np

def dvect(pos_0, pos_1, box, pbc):
    """Computes the shortest distance between pos_0 and pos_1 using box and pbc info.
    
    Keyword Arguments:
    pos_0, pos_1 -- vector positions or arrays of vector positions. The sizes and shapes of pos_0 and pos_1 have to be compatible with numpy broadcasting.
    box -- A Box instance.
    pbc -- Three boolean values indicating which box directions are periodic.     
    """
    
    #convert positions to np.arrays
    pos_0 = np.asarray(pos_0)
    pos_1 = np.asarray(pos_1)
    
    #get box values
    avect = box.avect
    bvect = box.bvect
    cvect = box.cvect 
    
    #compute the non-periodic distances between pos_0 and pos_1
    delta = pos_1 - pos_0
    if delta.ndim == 1:
        delta = delta[np.newaxis]
    
    #create iterators based on pbc
    check = [xrange(1), xrange(1), xrange(1)]
    for i in xrange(3):
        if pbc[i]:
            check[i] = xrange(-1, 2)    
    
    #Add all combinations of system vectors to delta to identify shortest d vector(s)
    d = delta.copy()
    for x in check[0]:
        for y in check[1]:
            for z in check[2]:
                test = delta + (x * avect + y * bvect + z * cvect)
                d = np.where(d.T[0]**2+d.T[1]**2+d.T[2]**2 <= test.T[0]**2+test.T[1]**2+test.T[2]**2, d.T, test.T).T
    
    if len(d) == 1:
        return d[0]
    else:
        return d