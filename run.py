import main, sys
from setup_problem import setup_problem
import time
from mpi4py import MPI
import numpy as np

start = time.time()

comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size

# grid spacing
nz,nr = 40,40

# tolerances: (default: 1e-5, 1e-50, 1e4, 1e4)
rtol    = 1e-6
atol    = 1e-10
runtime = 10
vl      = -200

for i in range(1, len(sys.argv)/2+1):
        key = 2*i-1
        val = 2*i
        if sys.argv[key] == '-nz':
            nz = int(sys.argv[val])
        elif sys.argv[key] == '-nr':
            nr = int(sys.argv[val])
        elif sys.argv[key] == '-rtol':
            rtol = float(sys.argv[val])
        elif sys.argv[key] == '-atol':
            atol = float(sys.argv[val])
        elif sys.argv[key] == '-rtime':
            runtime = float(sys.argv[val])
        elif sys.argv[key] == '-vl':
            vl = float(sys.argv[val])

if rank == 0:
    print('Mesh size is {} by {} = {} nodes'.format(nz,nr, nz*nr))
    print('Running on {} processors'.format(size))

if nr%size != 0:
    if rank == 0:
        print 'Error: Number of nodes in R must be divisible by number of tasks'
    sys.exit()


# setup problem
z,r_loc = setup_problem(nz,nr,vl,rank,size)

main.mod.initialize(rtol,atol)
#main.mod.view()

nsave     = 100
t = np.logspace(-3,np.log10(runtime),nsave)
main.mod.ne_save = np.zeros([len(z),len(r_loc),nsave],order='F')
main.mod.phi_save = np.zeros([len(z),len(r_loc),nsave],order='F')

main.mod.run(runtime,t)

if rank == 0:
    print('Elapsed time = {:.2f} sec'.format(time.time()-start))


if rank == 0:
    Istart = 0
else:
    Istart = 1

if rank == size-1:
    Iend = len(r_loc)
else:
    Iend = -1

r = np.zeros(nr,dtype='d')
comm.Allgather([r_loc[Istart:Iend],MPI.DOUBLE], [r, MPI.DOUBLE])

phi_loc = main.mod.phi
phi = np.zeros([nr,nz],dtype='d')
comm.Allgather([phi_loc[:,Istart:Iend],MPI.DOUBLE], [phi, MPI.DOUBLE])
phi = phi.T

ne_loc = main.mod.ne_pl
ne = np.zeros([nr,nz],dtype='d')
comm.Allgather([ne_loc[:,Istart:Iend],MPI.DOUBLE], [ne, MPI.DOUBLE])
ne = ne.T

ne_zt = main.mod.ne_save[:,0,:]
phi_zt = main.mod.phi_save[:,0,:]

ne_rt = main.mod.ne_save[0,:,:]
phi_rt = main.mod.phi_save[0,:,:]

if rank == 0:
    from plotter import *
    import matplotlib.pyplot as plt
    
    if nr > 1 and nz > 1:
        plot_phi(z,r,t,phi,phi_zt)
        plot_n(z,r,t,ne,ne_zt,'Electron')
        plt.show()
    
    elif nz > 1:
        plot_1dphi(z,t,phi_zt.T,'z')
        plot_1dn(z,t,ne_zt.T,'z','Electron')
        plt.show()
        
    elif nz > 1:
        plot_1dphi(z,t,phi_rt.T,'r')
        plot_1dn(z,t,ne_rt.T,'r','Electron')
        plt.show()
