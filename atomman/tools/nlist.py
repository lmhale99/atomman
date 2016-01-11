import numpy as np
from copy import deepcopy
from dvect import dvect
import time

def nlist(system, cutoff, cmult=1):
    start = time.time()
    natoms = system.natoms()
    vects, unit = system.box('vects')
    origin = system.box('origin')[0]
    pos = system.atoms.view['pos']
    pbc = system.pbc
    assert unit == system.atoms.unit['pos'],  'units do not match'
    if isinstance(cutoff, tuple) and len(cutoff) == 2:
        assert unit == cutoff[1],              'units do not match'
        cutoff = cutoff[0]
        
    #Determine orthogonal superbox that fully encompases the system
    corners = origin + np.array([[0,0,0], vects[0], vects[1], vects[2],
                                 vects[0] + vects[1],
                                 vects[0] + vects[2],
                                 vects[1] + vects[2],
                                 vects[0] + vects[1] + vects[2]])
    supermin = corners.min(axis=0) - cutoff
    supermax = corners.max(axis=0) + cutoff

    #Construct bins
    binsize = cutoff/cmult
    xbins = np.arange(supermin[0], supermax[0], binsize)
    ybins = np.arange(supermin[1], supermax[1], binsize)
    zbins = np.arange(supermin[2], supermax[2], binsize)
  
    #Create index lists for real atoms
    x_index = np.digitize(pos[:, 0], xbins) - 1
    y_index = np.digitize(pos[:, 1], ybins) - 1
    z_index = np.digitize(pos[:, 2], zbins) - 1
    xyz_index = np.hstack((x_index[:, np.newaxis], y_index[:, np.newaxis], z_index[:, np.newaxis]))
    atom_index = np.arange(natoms, dtype=int)
    
    real_bins = unique_rows(xyz_index)
    
    #create iterators based on pbc
    check = [xrange(1), xrange(1), xrange(1)]
    for i in xrange(3):
        if pbc[i]:
            check[i] = xrange(-1, 2) 

    #construct list of ghost atoms in the superbox
    
    for x in check[0]:
        for y in check[1]:
            for z in check[2]:
                if x == 0 and y == 0 and z == 0:
                    pass
                else:
                    newpos = x*vects[0] + y*vects[1] + z*vects[2] + pos
                    bool_check = np.empty_like(newpos, dtype=bool)
                    for i in xrange(3):
                        bool_check[:, i] = np.all([newpos[:, i] > supermin[i], newpos[:, i] < supermax[i]], axis=0)

                    new_index = np.where(np.all(bool_check, axis=1))[0]
                    try:
                        ghost_pos = np.vstack((ghost_pos, newpos[new_index]))
                        ghost_index = np.hstack((ghost_index, new_index))
                    except:
                        ghost_pos = newpos[new_index]
                        ghost_index = new_index

    #append index lists with ghost atoms
    if len(ghost_pos) > 0:
        x_index = np.digitize(ghost_pos[:, 0], xbins) - 1
        y_index = np.digitize(ghost_pos[:, 1], ybins) - 1
        z_index = np.digitize(ghost_pos[:, 2], zbins) - 1
        xyz_g_index = np.hstack((x_index[:, np.newaxis], y_index[:, np.newaxis], z_index[:, np.newaxis]))
        xyz_index = np.vstack((xyz_index, xyz_g_index))
        atom_index = np.hstack((atom_index, ghost_index))
    
    neighbors = np.zeros((natoms, 41), dtype=np.int)
    
    xyz_bins = np.zeros((len(xbins), len(ybins), len(zbins), 41), dtype=np.int)
    for i in xrange(len(atom_index)):
        x,y,z = xyz_index[i]
        xyz_bins[x,y,z,0] += 1
        try:
            xyz_bins[x,y,z, xyz_bins[x,y,z,0]] = atom_index[i]
        except:
            old_size = len(xyz_bins[0,0,0])
            newbins = np.zeros((len(xbins), len(ybins), len(zbins), old_size + 10), dtype=np.int)
            newbins[:, :, :, :old_size] = xyz_bins[:, :, :, :old_size]
            xyz_bins = newbins
            xyz_bins[x,y,z, xyz_bins[x,y,z,0]] = atom_index[i]
    
    end = time.time()
    print 'setup time:', end-start

    start = time.time()

    for bin in real_bins:
        x,y,z = bin
        
        short = xyz_bins[x,y,z, 1:xyz_bins[x,y,z,0]+1]
        long = deepcopy(short)
        
        for dx, dy, dz in box_iter(cmult):
            try:
                long = np.hstack((long, xyz_bins[x+dx,y+dy,z+dz, 1:xyz_bins[x+dx,y+dy,z+dz,0]+1]))
            except:
                long = xyz_bins[x+dx,y+dy,z+dz, 1:xyz_bins[x+dx,y+dy,z+dz,0]+1]
        
        for u in xrange(len(short)):
            u_pos = pos[short[u]]
            u_index = short[u]
            medium = long[u+1:]
            v_pos = pos[medium]
            try:
                d = np.linalg.norm(dvect(u_pos, v_pos, system.box(), system.pbc), axis=1)
            except:
                d = np.linalg.norm([dvect(u_pos, v_pos, system.box(), system.pbc)], axis=1)
            vlist = medium[np.where(d < cutoff)]
            for v_index in vlist:
                neighbors = append_neighbor(neighbors, u_index, v_index)
    
    for i in xrange(len(neighbors)):
        neighbors[i][1 : neighbors[i][0]+1] = np.sort(neighbors[i][1 : neighbors[i][0]+1])
    
    end = time.time()
    print 'neighbor assignment time:', end-start
    
    return neighbors
    
def append_neighbor(n, a, b):
    if b not in n[a, 1:n[a, 0]+1] and a != b:
        n[a, 0] += 1
        n[b, 0] += 1
        try:
            n[a, n[a, 0]] = b
            n[b, n[b, 0]] = a
        except:
            old_size = len(n[0])
            newlist = np.zeros((len(n), old_size + 10), dtype=np.int)
            newlist[:, :old_size] = n[:, :old_size]
            n = newlist
            n[a, n[a, 0]] = b
            n[b, n[b, 0]] = a
    return n
  
def unique_rows(a):  
    return np.unique(a.view(np.dtype((np.void, a.dtype.itemsize*a.shape[1])))).view(a.dtype).reshape(-1, a.shape[1])  
  
  
def box_iter(num):
    #Iterates over unique x,y,z combinations where each x,y,z can equal +-num.
    #Excludes 0,0,0 and equal but opposite combinations (i.e. x,y,z used, -x,-y,-z excluded)
    z = -num
    end = False
    while z <= 0:
        y = -num
        while y <= num:
            x = -num
            while x <= num:
                if x ==0 and y==0 and z==0:
                    end=True
                    break
                yield x,y,z
                x+=1
            if end:
                break
            y+=1
        if end:
            break
        z+=1     
        
def write_nlist(fname, n_list):
    #Write a neighbor list to a file
    
    with open(fname, 'w') as f:
        f.write('#Generated neighbor list\n')
        f.write('#index n_index_1 n_index_2 ...\n')
        for i in xrange(len(n_list)):
            f.write(str(i))
            for j in xrange(1, n_list[i, 0]+1):
                f.write(' ' + str(n_list[i, j]))
            f.write('\n')

def read_nlist(fname):
    #Read in a neighbor list from a file
    
    with open(fname, 'r') as f:
        natoms = 0
        max_n = 1
        
        #First pass determines number of atoms and max number of neighbors
        for line in f:
            terms = line.split()
            n_n = len(terms)
            if terms[0][0] != '#' and n_n > 0:
                natoms += 1
                if n_n > max_n: max_n = n_n
                
        n_list = np.zeros((natoms, max_n), dtype=np.int)
    
        #Second pass gets values
        for line in f:
            terms = line.split()
            if len(terms) > 0 and terms[0][0] != '#':
                i = int(terms[0])
                n_list[i, 0] = len(terms) - 1
                for j in xrange(1, len(terms)):
                    n_list[i, j] = terms[j]
    
    return n_list