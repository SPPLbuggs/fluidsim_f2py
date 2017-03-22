import numpy as np
import main
import sys

def setup_problem(nz,nr,v_l,rank,size):
    
    # Physical Params:
    length   = 5e-3
    radius   = 5e-3
    e_start  = 0.
    e_stop   = 5e-4
    
    # RF Params:
    v_r = 0
    rf_fq   = 1e6
    
    # Simulation Params:
    init = 1e14
    Te_init = 0.5
    len_0 = 1e-3
    phi_0 = 1e3
    tau_0 = 1e-6
    n_0 = 1e16
    
    if nz > 1 and nr > 1:
        stencil = 5
    else:
        stencil = 3
    
    # hyperbolic tangent spacing
    if nz > 1:
        z = np.linspace(-1.5,1.5,nz)
        z = np.tanh(z)*0.5+0.5
        z = z-z[0]
        z = z/z[-1]
        z = z*length/len_0
    else:
        z = np.array([0])
    
    if nr > 1:
        r = np.linspace(-1.5,0,nr)
        r = np.tanh(r)*0.5+0.5
        r = r-r[0]
        r = r/r[-1]
        r = r*radius/len_0
    else:
        r = np.array([0])
    
    # linear spacing
    #z = np.linspace(0,length,nz)/len_0
    #r = np.linspace(0,radius,nr)/len_0

    phi = np.zeros([nz,nr],order='F')/phi_0

    type_z = np.zeros([nz,nr],'int')
    type_r = np.zeros([nz,nr],'int')
    
    if nz > 1:
        for i in range(nz):
            for j in range(nr):
                if j == nr - 1:
                    type_r[i,j] = 1
                elif j == 0:
                    type_r[i,j] = -1
            
                if i == nz-1:
                    if e_start/len_0 <= r[j] <= e_stop/len_0:
                        type_z[i,j] = 2
                        phi[i,j] = v_r/phi_0
                    else:
                        type_z[i,j] = 1
                elif i == 0:
                    if e_start/len_0 <= r[j] <= e_stop/len_0:
                        type_z[i,j] = -2
                        phi[i,j] = v_l/phi_0
                    else:
                        type_z[i,j] = -1
    else:
        i = 0
        for j in range(nr):
            if j == 0:
                type_r[i,j] = -2
                phi[i,j] = v_l/phi_0
            if j == nr-1:
                type_r[i,j] =  2
                phi[i,j] = v_r/phi_0
    
    node_global = np.zeros([nz,nr],'int',order='F')
    
    neqn = 0
    for j in range(nr):
        for i in range(nz):
            if abs(type_z[i,j]) != 2 and abs(type_r[i,j]) != 2:
                neqn += 1
                node_global[i,j] = neqn
    
    Jstart = max(0,nr/size*rank-1)
    Jend   = min(nr,nr/size*(rank+1))+1
    
    
    r   = r[Jstart:Jend]
    phi = phi[:,Jstart:Jend]
    type_z = type_z[:,Jstart:Jend]
    type_r = type_r[:,Jstart:Jend]
    node_global = node_global[:,Jstart:Jend]
     
    ne = np.ones([nz,len(r)])*init/n_0
    ni = ne
    nt = ne*Te_init
    
    main.mod.neqn      = neqn
    main.mod.nz        = nz
    main.mod.nr        = nr
    main.mod.nr_loc    = len(r)
    main.mod.phir      = v_r
    main.mod.phil      = v_l
    
    main.mod.type_z    = type_z
    main.mod.type_r    = type_r
    main.mod.glob_node = node_global
    main.mod.z         = z
    main.mod.r         = r
    main.mod.phi       = phi
    main.mod.ne_pl     = ne
    main.mod.ni_pl     = ni
    main.mod.nt_pl     = nt
    
    return z,r
