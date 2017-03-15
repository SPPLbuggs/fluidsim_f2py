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

for i in range(1, len(sys.argv)/2+1):
        key = 2*i-1
        val = 2*i
        if sys.argv[key] == '-nz':
            nz = int(sys.argv[val])
        elif sys.argv[key] == '-nr':
            nr = int(sys.argv[val])

if rank == 0:
    print('Mesh size = {} by {}'.format(nz,nr))

if nz%size != 0:
    if rank == 0:
        print 'Error: Number of nodes in Z must be divisible by number of tasks'
    sys.exit()


# setup problem
r,z,phi = setup_problem(nz,nr,rank,size)

main.mod.initialize(1e-7,1e-7,5e5)
#main.mod.view()


runtime = 1
main.mod.run(runtime)

if rank == 0:
    print('Elapsed time = {:.2f} sec'.format(time.time()-start))


if rank == 0:
    Istart = 0
else:
    Istart = 1

if rank == size-1:
    Iend = len(z)
else:
    Iend = -1


z_global   = np.zeros(nz,dtype='d')
comm.Allgather([z[Istart:Iend].T,MPI.DOUBLE], [z_global, MPI.DOUBLE])

phi = main.mod.phi
phi_global = np.zeros([nz,nr],dtype='d')
comm.Allgather([phi[:,Istart:Iend].T,MPI.DOUBLE], [phi_global, MPI.DOUBLE])

ne = main.mod.ne_pl
ne_global = np.zeros([nz,nr],dtype='d')
comm.Allgather([ne[:,Istart:Iend].T,MPI.DOUBLE], [ne_global, MPI.DOUBLE])

if rank == 0:
    import matplotlib.pyplot as plt
    plt.figure()
    plt.contourf(r,z_global,phi_global,30)
    plt.colorbar()
    plt.xlabel('Radius (mm)')
    plt.ylabel('Z-axis (mm)')
    plt.title('Potential (V)')
    plt.show()
    
    plt.figure()
    plt.contourf(r,z_global,ne_global,30)
    plt.colorbar()
    plt.xlabel('Radius (mm)')
    plt.ylabel('Z-axis (mm)')
    plt.title('Electron Density ()')
    plt.show()

