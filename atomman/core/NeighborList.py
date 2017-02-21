#Standard library imports
from copy import deepcopy

#External library imports
import numpy as np

#Internal imports
from . import dvect
from . import System
from ..tools import uber_open_rmode

class NeighborList(object):
    """Class that finds and stores the neighbor atoms for a system."""
    
    def __init__(self, *args, **kwargs):
        """
        Class initializer. Either calls NeighborList.build() or NeighborList.load().
        
        Arguments:
        system -- atomman.System to calculate the neighbor list for.
        cutoff -- radial cutoff distance for neighbors.
        
        Keyword Arguments:
        cmult -- parameter associated with the binning routine.  Default value is most likely the fastest.
        initialsize -- specifies how many open neighbor positions to initially assign to each atom. Default value is 20.
        """
        if isinstance(args[0], System):
            self.build(*args, **kwargs)
        else:
            self.load(*args, **kwargs)
        
    
    @property
    def coord(self):
        """The atomic coordination numbers"""
        return deepcopy(self.__coord)
    
    def __len__(self):
        return len(self.__coord)
    
    def __getitem__(self, key):
        """get returns the list of neighbors for the specified atom."""
        if isinstance(key, (int, long, np.int)):
            return self.__neighbors[key, :self.coord[key]]
        else:
            raise IndexError('Index must be an int')
    
    def build(self, system, cutoff, cmult=1, initialsize=20):
        """
        Builds the neighbor list for a system.
        
        Arguments:
        system -- atomman.System to calculate the neighbor list for.
        cutoff -- radial cutoff distance for neighbors.
        
        Keyword Arguments:
        cmult -- parameter associated with the binning routine.  Default value is most likely the fastest.
        initialsize -- specifies how many open neighbor positions to initially assign to each atom. Default value is 20.
        """
        
        #Retrieve information from the system
        natoms = system.natoms
        vects = system.box.vects
        origin = system.box.origin
        pos = system.atoms.view['pos']
        pbc = system.pbc
        
        #Initialize class attributes
        self.__coord = np.zeros(natoms, dtype=int)
        self.__neighbors = np.empty((natoms, initialsize), dtype=int)
            
        #Determine orthogonal superbox that fully encompases the system
        corners = origin + np.array([[0,0,0], vects[0], vects[1], vects[2],
                                     vects[0] + vects[1],
                                     vects[0] + vects[2],
                                     vects[1] + vects[2],
                                     vects[0] + vects[1] + vects[2]])
        supermin = corners.min(axis=0) - 1.01 * cutoff
        supermax = corners.max(axis=0) + 1.01 * cutoff

        #Construct bins
        binsize = cutoff/cmult
        xbins = np.arange(supermin[0], supermax[0], binsize)
        ybins = np.arange(supermin[1], supermax[1], binsize)
        zbins = np.arange(supermin[2], supermax[2], binsize)

        #Build xyz box index for each atom
        x_index = np.digitize(pos[:, 0], xbins) - 1
        y_index = np.digitize(pos[:, 1], ybins) - 1
        z_index = np.digitize(pos[:, 2], zbins) - 1
        xyz_index = np.hstack((x_index[:, np.newaxis], y_index[:, np.newaxis], z_index[:, np.newaxis]))
        
        #Relate atom's id to xyz_index
        atom_index = np.arange(natoms, dtype=int)
        
        #Identify all bins with real atoms
        real_bins = self.__unique_rows(xyz_index)
        
        #create iterators based on pbc
        check = [(0), (0), (0)]
        for i in xrange(3):
            if pbc[i]:
                check[i] = (-1, 0, 1) 

        #construct list of ghost atoms in the superbox
        ghost_pos = []
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
        
        
        #assign atoms and ghost atoms to xyz bins
        xyz_bins = np.zeros((len(xbins), len(ybins), len(zbins), initialsize*2), dtype=np.int)
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
        
        #iterate over all bins with real atoms
        for bin in real_bins:
            x,y,z = bin
            
            #short = all atoms in current bin
            short = xyz_bins[x,y,z, 1:xyz_bins[x,y,z,0]+1]
            if x == 1 and y == 1 and z == 0:
                print short
            long = deepcopy(short)
            
            #add all atoms in half of the nearby bins to long
            for dx, dy, dz in self.__box_iter(cmult):
                try:
                    long = np.hstack((long, xyz_bins[x+dx,y+dy,z+dz, 1:xyz_bins[x+dx,y+dy,z+dz,0]+1]))
                except:
                    long = xyz_bins[x+dx,y+dy,z+dz, 1:xyz_bins[x+dx,y+dy,z+dz,0]+1]
            
            #compare all atoms in short to medium (which is long starting with u+1).
            for u in xrange(len(short)):
                u_pos = pos[short[u]]
                u_index = short[u]
                medium = long[u+1:]
                v_pos = pos[medium]
                try:
                    d = np.linalg.norm(dvect(u_pos, v_pos, system.box, system.pbc), axis=1)
                except:
                    d = np.linalg.norm([dvect(u_pos, v_pos, system.box, system.pbc)], axis=1)
                vlist = medium[np.where(d < cutoff)]
                for v_index in vlist:
                    self.__append_neighbor(u_index, v_index)
        
        #sort each atom's neighbors
        for i in xrange(len(self.__neighbors)):
            self.__neighbors[i, :self.coord[i]] = np.sort(self.__neighbors[i, :self.coord[i]])
    
    def __append_neighbor(self, a, b):
        """Adds atom ids a and b to each other's list of neighbors."""
        
        if b not in self[a] and a != b:
            
            #Increase coordination number
            self.__coord[a] += 1
            self.__coord[b] += 1
            
            #Add a and b to each other's lists
            try:
                self.__neighbors[a, self.coord[a]-1] = b
                self.__neighbors[b, self.coord[b]-1] = a
            
            #Increase list size if needed
            except:
                oldsize = len(self.__neighbors[0])
                newlist = np.empty((len(self.__neighbors), oldsize + 5), dtype=int)
                newlist[:, :oldsize] = self.__neighbors[:, :oldsize]
                self.__neighbors = newlist
                self.__neighbors[a, self.coord[a]-1] = b
                self.__neighbors[b, self.coord[b]-1] = a
  
    def __unique_rows(self, a):  
        """Takes two-dimensional array a and returns only the unique rows."""
        return np.unique(a.view(np.dtype((np.void, a.dtype.itemsize*a.shape[1])))).view(a.dtype).reshape(-1, a.shape[1])  
  
    def __box_iter(self, num):
        """
        Iterates over unique x,y,z combinations where each x,y,z can equal +-num.
        Excludes 0,0,0 and equal but opposite combinations (e.g. x,y,z used, -x,-y,-z excluded).
        """
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
            
    def load(self, fp, initialsize=20):
        """Read in a neighbor list from a file"""
        
        #First pass determines number of atoms and max number of neighbors
        with uber_open_rmode(fp) as fin:
            for line in fin:
                terms = line.split()
                n_n = len(terms)
                if terms[0][0] != '#' and n_n > 0:
                    natoms += 1
                    if n_n > initialsize: initialsize = n_n
            
            self.__coord = np.zeros(natoms, dtype=int)
            self.__neighbors = np.empty((natoms, initialsize), dtype=int)

            #Second pass gets values
            for line in fin:
                terms = line.split()
                if len(terms) > 0 and terms[0][0] != '#':
                    i = int(terms[0])
                    self.__coord[i] = len(terms) - 1
                    for j in xrange(1, len(terms)):
                        self.__neighbors[i, j-1] = terms[j]
    
    def dump(self, fname):
        """Saves the neighbor list to a file"""
        with open(fname, 'w') as fp:
            fp.write('# Neighbor list:\n')
            fp.write('# The first column gives an atom index.\n')
            fp.write('# The rest of the columns are the indexes of the identified neighbors.\n')
            for i in xrange(len(self)):
                fp.write('%i' % i)
                for j in self[i]:
                    fp.write(' %i' % j)
                fp.write('\n')
            
        
    