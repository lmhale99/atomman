import numpy as np
import atomman as am

def nye_tensor(system, p_vectors, theta_max = 27, axes=None, neighbor_list=None, neighbor_list_cutoff=None):
    """Computes strain properties and Nye tensor for a defect containing system."""

    #neighbor list setup
    if neighbor_list is not None:
        assert neighbor_list_cutoff is None, 'neighbor_list and neighbor_list_cutoff cannot both be given'
    elif neighbor_list_cutoff is not None:
        neighbor_list = am.nlist(system, neighbor_list_cutoff)
    elif 'nlist' in system.prop:
        neighbor_list = system.prop['nlist']
    
    #If p_vectors is only given for one atom, apply to all atoms
    if p_vectors.ndim == 2:
        p_vectors = np.broadcast_to(p_vectors, (system.natoms, len(p_vectors), 3))
    
    #If axes are given, transform p accordingly
    if axes is not None:
        p_vectors = np.inner(p_vectors, am.tools.axes_check(axes))
        
    #get cos of theta_max
    cos_theta_max = np.cos(theta_max * np.pi / 180)  
    
    #Define identity and epsilon arrays   
    #iden = 
    eps = np.array([[[ 0, 0, 0],[ 0, 0, 1],[ 0,-1, 0]],
                    [[ 0, 0,-1],[ 0, 0, 0],[ 1, 0, 0]],
                    [[ 0, 1, 0],[-1, 0, 0],[ 0, 0, 0]]])
    
    #Identify largest number of nearest neighbors
    nmax = 0
    for ns in neighbor_list:
        if ns[0] > nmax: nmax=ns[0]
    
    #Initialize variables (done here to reduce memory allocations making it slightly faster)
    strain = np.empty((system.natoms, 3, 3))
    inv1 = np.empty(system.natoms)
    inv2 = np.empty(system.natoms)
    ang_vel = np.empty(system.natoms)
    nye = np.empty((system.natoms, 3, 3))
    P = np.zeros((nmax, 3))            
    Q = np.zeros((nmax, 3))
    G = np.empty((system.natoms, 3, 3))
    gradG = np.empty((3, 3, 3))
    
    #Loop to calculate correspondence tensor, G, and strain data   
    for i in xrange(system.natoms):
        p = p_vectors[i]
        p_mags = np.linalg.norm(p, axis=1)
        r1 = p_mags.min()
        
        #Calculate radial neighbor vectors, q
        q = system.dvect(i, neighbor_list[i][1:neighbor_list[i][0]+1])
        q_mags = np.linalg.norm(q, axis=1)
        
        #Calculate cos_thetas between all p's and q's.
        cos_thetas = (np.dot(p, q.T) /q_mags ).T / p_mags
        
        #Identify best index matches 
        index_pairing = cos_thetas.argmax(1)
        
        #Exclude values where theta is greater than theta_max
        index_pairing[cos_thetas.max(1) < cos_theta_max] = -1
        
        #Search for duplicate index_pairings
        u, u_count = np.unique(index_pairing, return_counts=True)
        for match in u[(u!=-1) & (u_count > 1)]:
            print index_pairing==match
        
        for n in xrange(neighbor_list[i][0]):               
            #Check if the particular p has already been assigned to another q
            #Remove the p-q pair that is farther from r1
            if index_pairing[n] >=0:
                for k in xrange(n):
                    if index_pairing[n] == index_pairing[k]:
                        nrad = abs(r1 - q_mags[n])
                        krad = abs(r1 - q_mags[k])
                        if nrad < krad:
                            index_pairing[k]=-1
                        else:
                            index_pairing[n]=-1

        #Construct reduced P, Q matrices from p-q pairs
        c = 0
        for n in xrange(neighbor_list[i][0]):
            if index_pairing[n] >= 0:
                Q[c] = q[n]
                P[c] = p[index_pairing[n]]
                c+=1   
        #Compute lattice correspondence tensor, G, from P and Q
        if c == 0:
            G[i] = np.identity(3)
            print i+1
            raise ValueError('An atom lacks pair sets. Check neighbor list')
        else:   
            G[i], resid, rank, s = np.linalg.lstsq(Q[:c], P[:c])

        #Compute strain properties
   
        strain[i] = ((np.identity(3) - G[i]) + (np.identity(3) - G[i]).T) / 2.
        inv1[i] = strain[i,0,0] + strain[i,1,1] + strain[i,2,2]
        inv2[i] = (strain[i,0,0] * strain[i,1,1] + strain[i,0,0] * strain[i,2,2] + strain[i,1,1] * strain[i,2,2] 
                - strain[i,0,1]**2 - strain[i,0,2]**2 - strain[i,1,2]**2)
        rot = ((np.identity(3) - G[i]) - (np.identity(3) - G[i]).T) / 2.
        ang_vel[i] = (rot[0,1]**2 + rot[0,2]**2 + rot[1,2]**2)**0.5

    #Construct the gradient tensor of G, gradG
    for i in xrange(system.natoms):
        neighbors = neighbor_list[i][1:neighbor_list[i][0]+1]
        Q = system.dvect(i, neighbors)
        dG = G[neighbors] - G[i]        
        for x in xrange(3):
            gradG[x,:] = np.linalg.lstsq(Q, dG[:,x,:])[0].T
        
        #Use gradG to build the nye tensor, nye
        nye[i] = -np.einsum('ijm,ikm->jk', eps, gradG) 
    
    return {'strain':strain, 'strain_invariant_1':inv1, 'strain_invariant_2':inv2, 'angular_velocity':ang_vel, 'Nye_tensor':nye}
