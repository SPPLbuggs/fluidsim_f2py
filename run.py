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
simtype = 'dc'

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
        elif sys.argv[key] == '-type':
            simtype = sys.argv[val]
        elif sys.argv[key] == '-vl':
            vl = float(sys.argv[val])

if simtype.lower() == 'dc':
    wvfm = 0
elif simtype.lower() == 'ac':
    wvfm = 1
elif simtype.lower() == 'pulse':
    wvfm = 2
else:
    if rank == 0:
        print 'Error: unknown simulation type. Options are dc, ac, and pulse.'
    sys.exit()

if rank == 0:
    print('Mesh size is {} by {} = {} nodes'.format(nz,nr, nz*nr))
    print('Running on {} processors'.format(size))

if nr%size != 0:
    if rank == 0:
        print 'Error: Number of nodes in R must be divisible by number of tasks'
    sys.exit()


# setup problem
z,r_loc = setup_problem(nz,nr,vl,wvfm,rank,size)

main.mod.initialize(rtol,atol)
#main.mod.view()

nsave = 1000
#t = np.logspace(-3,np.log10(runtime),nsave)
t = np.linspace(1,runtime,nsave)
main.mod.phi_save = np.zeros([len(z),len(r_loc),nsave],order='F')
main.mod.ni_save = np.zeros([len(z),len(r_loc),nsave],order='F')
main.mod.ne_save = np.zeros([len(z),len(r_loc),nsave],order='F')
main.mod.nt_save = np.zeros([len(z),len(r_loc),nsave],order='F')
main.mod.nm_save = np.zeros([len(z),len(r_loc),nsave],order='F')

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

ni_loc = main.mod.ni_pl
ni = np.zeros([nr,nz],dtype='d')
comm.Allgather([ni_loc[:,Istart:Iend],MPI.DOUBLE], [ni, MPI.DOUBLE])
ni = ni.T

ne_loc = main.mod.ne_pl
ne = np.zeros([nr,nz],dtype='d')
comm.Allgather([ne_loc[:,Istart:Iend],MPI.DOUBLE], [ne, MPI.DOUBLE])
ne = ne.T

nt_loc = main.mod.nt_pl
nt = np.zeros([nr,nz],dtype='d')
comm.Allgather([nt_loc[:,Istart:Iend],MPI.DOUBLE], [nt, MPI.DOUBLE])
nt = nt.T

nm_loc = main.mod.nm_pl
nm = np.zeros([nr,nz],dtype='d')
comm.Allgather([nm_loc[:,Istart:Iend],MPI.DOUBLE], [nm, MPI.DOUBLE])
nm = nm.T

phi_zt = main.mod.phi_save[:,0,:]
ni_zt = main.mod.ni_save[:,0,:]
ne_zt = main.mod.ne_save[:,0,:]
nt_zt = main.mod.nt_save[:,0,:]
nm_zt = main.mod.nm_save[:,0,:]

phi_rt = main.mod.phi_save[0,:,:]
ni_rt = main.mod.ni_save[0,:,:]
ne_rt = main.mod.ne_save[0,:,:]
nt_rt = main.mod.nt_save[0,:,:]
nm_rt = main.mod.nm_save[0,:,:]

if rank == 0:
    from plotter import *
    import matplotlib.pyplot as plt
    
    if nr > 1 and nz > 1:
        plot_phi(z,r,t,phi,phi_zt)
        plot_n(z,r,t,ne,ne_zt,'Electron')
        plot_n(z,r,t,ni,ni_zt,'Ion')
        plot_n(z,r,t,nm,nm_zt,'Metastable')
        plot_t(z,r,t,nt/ne,nt_zt/ne_zt)
        plt.show()
    
    elif nz > 1:
        plot_1dphi(z,t,phi_zt.T,'z')
        plot_1dn(z,t,ne_zt.T,'z','Electron')
        plot_1dn(z,t,ni_zt.T,'z','Ion')
        plot_1dn(z,t,nm_zt.T,'z','Metastable')
        plot_1dt(z,t,nt_zt.T/ne_zt.T,'z')
        plt.show()
        
    elif nr > 1:
        plot_1dphi(r,t,phi_rt.T,'r')
        plot_1dn(r,t,ne_rt.T,'r','Electron')
        plot_1dn(r,t,ni_rt.T,'r','Ion')
        plot_1dn(r,t,nm_rt.T,'r','Metastable')
        plot_1dt(r,t,nt_rt.T/ne_zt.T,'r')
        plt.show()
