import numpy as np
from atomman.tools import mag, axes_check

def nye(system, p, neighbor_list=None, neighbor_list_cutoff=None, axes=None, tmax = 27):
    axes = system.prop('axes')
    nlist = system.prop('nlist')
    natoms = system.natoms()
    
    T,vmag = axes_check(axes)
    if len(p) == 1:
        p = p[0]
    else:
        raise ValueError('Multiple unique sites not supported yet')
    
    p = p.dot(np.transpose(T))
    
    #Change tmax=angle to tmax=cos(angle)
    tmax = np.cos(tmax * np.pi / 180)  
    iden = np.identity(3)
    eps = np.array([[[ 0, 0, 0],[ 0, 0, 1],[ 0,-1, 0]],
                    [[ 0, 0,-1],[ 0, 0, 0],[ 1, 0, 0]],
                    [[ 0, 1, 0],[-1, 0, 0],[ 0, 0, 0]]])
    
    #Set r1 to smallest interatomic distance from p
    r1 = 50.
    for pi in p:
        if r1 > mag(pi): r1 = mag(pi)
    
    #Identify largest number of nearest neighbors
    nmax = 0
    for ns in nlist:
        if ns[0] > nmax: nmax=ns[0]
    
    #Initialize variables (done here to reduce memory allocations making it slightly faster)
    q = np.zeros((nmax, 3))
    temp = [-1 for y in xrange(nmax)]
    P = np.zeros((nmax, 3))            
    Q = np.zeros((nmax, 3))
    G = [[] for y in xrange(natoms)]
    nye = np.zeros((3, 3))
    gradG = np.zeros((3, 3, 3))
    
    #Loop to calculate correspondence tensor, G, and strain data   
    for i in xrange(natoms):
        #Reset temp list used to identify p-q pairs
        for t in xrange(len(temp)):
            temp[t] = -1
        
        #Calculate radial neighbor vectors, q, and associate to perfect crystal vectors, p using theta.
        #Only match p-q if cos(angle) between them is largest and greater than tmax.
        for n in xrange(nlist[i][0]):
            j = nlist[i][n+1]
            q[n] = system.dvect(i, j)
            theta = tmax
            for s in xrange(len(p)):
                if (p[s].dot(q[n]) / (mag(p[s]) * mag(q[n])) > theta):
                    temp[n] = s
                    theta = p[s].dot(q[n]) / (mag(p[s]) * mag(q[n]))
                       
            #Check if the particular p has already been assigned to another q
            #Remove the p-q pair that is farther from r1
            if temp[n] >=0:
                for k in xrange(n):
                    if temp[n] == temp[k]:
                        nrad = abs(r1 - mag(q[n]))
                        krad = abs(r1 - mag(q[k]))
                        if nrad < krad:
                            temp[k]=-1
                        else:
                            temp[n]=-1

        #Construct reduced P, Q matrices from p-q pairs
        c = 0
        for n in xrange(nlist[i][0]):
            if temp[n] >= 0:
                Q[c] = q[n]
                P[c] = p[temp[n]]
                c+=1   
        #Compute lattice correspondence tensor, G, from P and Q
        if c == 0:
            G[i] = iden
            print i+1
            raise ValueError('An atom lacks pair sets. Check neighbor list')
        else:   
            G[i], resid, rank, s = np.linalg.lstsq(Q[:c], P[:c])
        
        #Compute strain properties
        strain = ((iden - G[i]) + (iden - G[i]).T) / 2.
        inv1 = strain[0,0] + strain[1,1] + strain[2,2]
        inv2 = (strain[0,0] * strain[1,1] + strain[0,0] * strain[2,2] + strain[1,1] * strain[2,2] 
                - strain[0,1]**2 - strain[0,2]**2 - strain[1,2]**2)
        rot = ((iden - G[i]) - (iden - G[i]).T) / 2.
        ang_vel = (rot[0,1]**2 + rot[0,2]**2 + rot[1,2]**2)**0.5
        
        system.atoms(i, 'strain', strain)
        system.atoms(i, 'strain_invariant_1', inv1)
        system.atoms(i, 'strain_invariant_2', inv2)
        system.atoms(i, 'angular_velocity', ang_vel)

    #Construct the gradient tensor of G, gradG
    for i in xrange(natoms):
        for x in xrange(3):
            for y in xrange(3):
                Q = np.zeros((nlist[i][0], 3))
                dG = np.zeros((nlist[i][0]))
                for n in xrange(nlist[i][0]):
                    j = nlist[i][n+1]
                    Q[n] = system.dvect(i, j)
                    dG[n] = G[j][x,y] - G[i][x,y]
                gradG[x, y], resid, rank, s = np.linalg.lstsq(Q, dG)  
        
        #Use gradG to build the nye tensor, nye
        for ij in xrange(3):
            nye[0,ij] = -gradG[2, ij, 1] + gradG[1, ij, 2]
            nye[1,ij] = -gradG[0, ij, 2] + gradG[2, ij, 0]
            nye[2,ij] = -gradG[1, ij, 0] + gradG[0, ij, 1]
                
        system.atoms(i, 'Nye', nye)