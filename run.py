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
rtol = 1e-6
atol = 1e-10
runtime = 10

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

if rank == 0:
    print('Mesh size is {} by {} = {} nodes'.format(nz,nr, nz*nr))
    print('Running on {} processors'.format(size))

if nr%size != 0:
    if rank == 0:
        print 'Error: Number of nodes in R must be divisible by number of tasks'
    sys.exit()


# setup problem
z,r,phi = setup_problem(nz,nr,rank,size)

main.mod.initialize(rtol,atol)
#main.mod.view()

nsave     = 100
savetimes = np.logspace(-3,np.log10(runtime),nsave)
main.mod.ne_save = np.zeros([len(z),len(r),nsave],order='F')
main.mod.phi_save = np.zeros([len(z),len(r),nsave],order='F')

main.mod.run(runtime,savetimes)

if rank == 0:
    print('Elapsed time = {:.2f} sec'.format(time.time()-start))


if rank == 0:
    Istart = 0
else:
    Istart = 1

if rank == size-1:
    Iend = len(r)
else:
    Iend = -1

r_global   = np.zeros(nr,dtype='d')
comm.Allgather([r[Istart:Iend],MPI.DOUBLE], [r_global, MPI.DOUBLE])

phi = main.mod.phi
phi_global = np.zeros([nr,nz],dtype='d')
comm.Allgather([phi[:,Istart:Iend],MPI.DOUBLE], [phi_global, MPI.DOUBLE])

ne = main.mod.ne_pl
ne_global = np.zeros([nr,nz],dtype='d')
comm.Allgather([ne[:,Istart:Iend],MPI.DOUBLE], [ne_global, MPI.DOUBLE])

ne_zt = main.mod.ne_save[:,0,:]
phi_zt = main.mod.phi_save[:,0,:]

if rank == 0:
    from plotter import *
    import matplotlib.pyplot as plt
    
    if nr > 1:
        plot_phi(z,r,savetimes,phi_global.T,phi_zt)
        plot_n(z,r,savetimes,ne_global.T,ne_zt,'Electron')
        plt.show()
    
    else:
        plot_1dphi(z,savetimes,phi_zt.T)
        plot_1dn(z,savetimes,ne_zt.T)
        plt.show()

