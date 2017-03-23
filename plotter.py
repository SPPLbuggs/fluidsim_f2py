import numpy as np
import matplotlib.pyplot as plt

size = 12
med_size = 13
big_size = 14

plt.rc('font', size=size)
plt.rc('axes', titlesize=size)
plt.rc('axes', labelsize=med_size)
plt.rc('xtick', labelsize=size)
plt.rc('ytick', labelsize=size)
plt.rc('legend', fontsize=size)
plt.rc('figure', titlesize=big_size)
plt.rcParams['figure.figsize'] = (4.5, 3)

def plot_phi(z,r,t,n,nt):

    rr,zz = np.meshgrid(r,z)
    tt,zt = np.meshgrid(t,z)
    
    fig,axes = plt.subplots(1,2,sharey=True,figsize=[8.3,3.2])
    (ax1,ax2) = axes
    vmin = nt.min()*1e3
    vmax = nt.max()*1e3
    ax1.contourf(rr,zz,n*1e3,30,vmin=vmin,vmax=vmax)
    ax1.contourf(-rr,zz,n*1e3,30,vmin=vmin,vmax=vmax)
    ax1.set_title('t = {:.2f} $\mu$s'.format(t[-1]))
    ax1.set_ylabel('z (mm)')
    ax1.set_xlabel('radius (mm)')
    im = ax2.contourf(tt,zt,nt*1e3,30,vmin=vmin,vmax=vmax)
    ax2.set_title('r = 0 mm')
    ax2.set_xlabel('time ($\mu$s)')
    plt.suptitle('Potential ($\phi$)')
    plt.subplots_adjust(left=0.08,right=1.05,top=0.85,bottom=0.2,wspace=0.125)
    clb = fig.colorbar(im, ax = axes.ravel().tolist())
    clb.ax.set_title('(V)')
    # plt.savefig('phi_{}V_{}.png'.format(V,RF),dpi=600)

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=[8,3])
    ax1.plot(z,n[:,len(r)/2]*1e3,label='r = {:.1f}'.format(z[len(z)/2]))
    ax1.plot(z,n[:,0]*1e3,label='r = {:.1f}'.format(z[0]))
    ax2.plot(r,n[-1,:]*1e3,label='anode')
    ax2.plot(r,n[0,:]*1e3,label='cathode')
    ax1.legend(frameon=False, loc='best')
    ax2.legend(frameon=False, loc='best')  
    ax1.set_xlabel('z (mm)')
    ax2.set_xlabel('radius (mm)')
    ax1.set_ylabel('Potential (V)')
    plt.subplots_adjust(left=0.12,right=0.96,bottom=0.2,top=0.9,wspace=0.125)

def plot_n(z,r,t,n,nt,particle):
    
    rr,zz = np.meshgrid(r,z)
    tt,zt = np.meshgrid(t,z)
    
    fig,axes = plt.subplots(1,2,sharey=True,figsize=[8.3,3.2])
    (ax1,ax2) = axes
    vmin = nt.min()
    vmax = nt.max()
    ax1.contourf(rr,zz,n,30,vmin=vmin,vmax=vmax)
    ax1.contourf(-rr,zz,n,30,vmin=vmin,vmax=vmax)
    ax1.set_title('t = {:.2f} $\mu$s'.format(t[-1]))
    ax1.set_ylabel('z (mm)')
    ax1.set_xlabel('radius (mm)')
    im = ax2.contourf(tt,zt,nt,30,vmin=vmin,vmax=vmax)
    ax2.set_title('r = 0 mm')
    ax2.set_xlabel('time ($\mu$s)')
    plt.suptitle(particle + ' Density')
    plt.subplots_adjust(left=0.08,right=1.05,top=0.85,bottom=0.2,wspace=0.125)
    clb = fig.colorbar(im, ax = axes.ravel().tolist())
    clb.ax.set_title('(10$^{15}$ m$^{-3}$)')
    # plt.savefig('phi_{}V_{}.png'.format(V,RF),dpi=600)

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=[8,3])
    ax1.plot(z,n[:,len(r)/2]*1e15,label='r = {:.1f}'.format(z[len(z)/2]))
    ax1.plot(z,n[:,0]*1e15,label='r = {:.1f}'.format(z[0]))
    ax2.plot(r,n[-1,:]*1e15,label='anode')
    ax2.plot(r,n[0,:]*1e15,label='cathode')
    ax1.legend(frameon=False, loc='best')
    ax2.legend(frameon=False, loc='best')  
    ax1.set_xlabel('z-axis (mm)')
    ax2.set_xlabel('r-axis (mm)')
    ax1.set_ylabel(particle + 'Density (m-3)')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    plt.subplots_adjust(left=0.12,right=0.96,top=0.9,bottom=0.2,wspace=0.125)

def plot_t(z,r,t,n,nt):
    
    rr,zz = np.meshgrid(r,z)
    tt,zt = np.meshgrid(t,z)
    
    fig,axes = plt.subplots(1,2,sharey=True,figsize=[8.3,3.2])
    (ax1,ax2) = axes
    vmin = nt.min()
    vmax = nt.max()
    ax1.contourf(rr,zz,n,30,vmin=vmin,vmax=vmax)
    ax1.contourf(-rr,zz,n,30,vmin=vmin,vmax=vmax)
    ax1.set_title('t = {:.2f} $\mu$s'.format(t[-1]))
    ax1.set_ylabel('z (mm)')
    ax1.set_xlabel('radius (mm)')
    im = ax2.contourf(tt,zt,nt,30,vmin=vmin,vmax=vmax)
    ax2.set_title('r = 0 mm')
    ax2.set_xlabel('time ($\mu$s)')
    plt.suptitle('Electron Temperature (T$_e$)')
    plt.subplots_adjust(left=0.08,right=1.05,top=0.85,bottom=0.2,wspace=0.125)
    clb = fig.colorbar(im, ax = axes.ravel().tolist())
    clb.ax.set_title('(eV)')
    # plt.savefig('phi_{}V_{}.png'.format(V,RF),dpi=600)

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=[8,3])
    ax1.plot(z,n[:,len(r)/2],label='r = {:.1f}'.format(z[len(z)/2]))
    ax1.plot(z,n[:,0],label='r = {:.1f}'.format(z[0]))
    ax2.plot(r,n[-1,:],label='anode')
    ax2.plot(r,n[0,:],label='cathode')
    ax1.legend(frameon=False, loc='best')
    ax2.legend(frameon=False, loc='best')  
    ax1.set_xlabel('z-axis (mm)')
    ax2.set_xlabel('r-axis (mm)')
    ax1.set_ylabel('Electron Temperature (eV)')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    plt.subplots_adjust(left=0.12,right=0.96,top=0.9,bottom=0.2,wspace=0.125)

def plot_1dphi(z,t,n,dir):
    zz,tt = np.meshgrid(z,t)
    
    fig,axes = plt.subplots(1,2,sharey=True,figsize=[8,3])
    (ax1,ax2) = axes

    ax1.plot(n[1,:]*1e3,z,label='{:.1f}us'.format(t[1]))
    ax1.plot(n[len(t)/2,:]*1e3,z,label='{:.1f}us'.format(t[len(t)/2]))
    ax1.plot(n[-1,:]*1e3,z,label='{:.1f}us'.format(t[-1]))
    ax1.set_xlabel('Potential (V)')
    ax1.set_ylabel(dir+' (mm)')
    ax1.legend(frameon=False, loc='best')
    
    im = ax2.contourf(tt,zz,n*1e3,30)
    ax2.set_xlabel('Time ($\mu$s)')
    
    plt.suptitle('Potential ($\phi$)')
    plt.subplots_adjust(left=0.08,right=1.05,top=0.85,wspace=0.125,bottom=0.2)
    clb = fig.colorbar(im, ax = axes.ravel().tolist())
    clb.ax.set_title('(V)')

def plot_1dn(z,t,n,dir,particle):
    zz,tt = np.meshgrid(z,t)
    
    fig,axes = plt.subplots(1,2,sharey=True,figsize=[8,3])
    (ax1,ax2) = axes

    ax1.plot(n[1,:]*1e15,z,label='{:.1f}us'.format(t[1]))
    ax1.plot(n[len(t)/2,:]*1e15,z,label='{:.1f}us'.format(t[len(t)/2]))
    ax1.plot(n[-1,:]*1e15,z,label='{:.1f}us'.format(t[-1]))
    ax1.set_xlabel(particle + ' Density (m-3)')
    ax1.set_ylabel(dir+' (mm)')
    ax1.set_xscale('log')
    ax1.legend(frameon=False, loc='best')
    
    im = ax2.contourf(tt,zz,n,30)
    ax2.set_xlabel('Time ($\mu$s)')
    
    plt.suptitle(particle + ' Density')
    plt.subplots_adjust(left=0.08,right=1.05,top=0.85,bottom=0.2,wspace=0.125)
    clb = fig.colorbar(im, ax = axes.ravel().tolist())
    clb.ax.set_title('(10$^{15}$ m$^{-3}$)')

def plot_1dt(z,t,n,dir):
    zz,tt = np.meshgrid(z,t)
    
    fig,axes = plt.subplots(1,2,sharey=True,figsize=[8,3])
    (ax1,ax2) = axes

    ax1.plot(n[1,:],z,label='{:.1f}us'.format(t[1]))
    ax1.plot(n[len(t)/2,:],z,label='{:.1f}us'.format(t[len(t)/2]))
    ax1.plot(n[-1,:],z,label='{:.1f}us'.format(t[-1]))
    ax1.set_xlabel('Electron Temperature (eV)')
    ax1.set_ylabel(dir+' (mm)')
    ax1.legend(frameon=False, loc='best')
    
    im = ax2.contourf(tt,zz,n,30)
    ax2.set_xlabel('Time ($\mu$s)')
    
    plt.suptitle('Electron Temperature (T$_e$)')
    plt.subplots_adjust(left=0.08,right=1.05,top=0.85,bottom=0.2,wspace=0.125)
    clb = fig.colorbar(im, ax = axes.ravel().tolist())
    clb.ax.set_title('(eV)')


        
