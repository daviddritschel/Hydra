#!/usr/bin/env python

# This script plots |u|^2 and |B|^2 at three selected times

#========== Perform the generic imports =========
import warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.it'] = 'STIXGeneral:italic'
mpl.rcParams['mathtext.bf'] = 'STIXGeneral:italic:bold'

# set tick label size:
label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
# set axes width:
mpl.rcParams['axes.linewidth'] = 2

#====================== Function definitions =======================
def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

    emag=1.0
    rmult=fmax-fmin
    while rmult < 10:
       emag=emag/10
       rmult=rmult*10

    while rmult >= 100:
       emag=emag*10
       rmult=rmult/10

    kmult=int(rmult/10)

    if kmult < 1:
       ci=emag
    elif kmult < 2:
       ci=2*emag
    elif kmult < 4:
       ci=4*emag
    elif kmult < 8:
       ci=10*emag
    else:
       ci=20*emag

    return ci

#==============================================================================
# Set up figure:
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(figsize=[20,13.3], nrows=2, ncols=3)

# Customise tick values:
ax1.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax4.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax5.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax6.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax4.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax5.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax6.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

ax1.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax4.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax5.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax6.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax4.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax5.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax6.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

#=================================================================
# Specify grid resolution:
print()
ng_in = input(' Grid resolution (default 512)? ')
ng = int(ng_in or 512)
N=ng*ng
print()

# Open ecomp.asc file in one directory to get time between frames:
in_file=open('evolution/ecomp.asc','r')
time, e1, e2, e3 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

dt=time[1]-time[0]
tmax=time[-1]

# Read data into arrays for plotting:
# -----------------------------------------------------------
t_def=0.1*tmax
t_in = input(' First time to show (default '+str(t_def)+')? ')
t1 = float(t_in or t_def)
frame=int((t1+0.0001)/dt)

Z=np.empty([ng+1,ng+1])
in_file=open('evolution/usq.r4','r')
raw_uarray=np.fromfile(in_file,dtype=np.float32)
in_file.close()

Z[0:ng,0:ng]=raw_uarray[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z[ng,0:ng]=Z[0,0:ng]
Z[0:ng+1,ng]=Z[0:ng+1,0]

zmax=np.amax(Z)
dz=contint(0.0,zmax)
jmin=0
jmax=int(zmax/dz)
clevels=np.linspace(0.0,dz*float(jmax),jmax-jmin+1)

im1=ax1.imshow(Z,cmap=cm.seismic,vmin=0.0,vmax=zmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im1, cax=cax, ticks=clevels)

Z=np.empty([ng+1,ng+1])
in_file=open('evolution/bsq.r4','r')
raw_barray=np.fromfile(in_file,dtype=np.float32)
in_file.close()

Z[0:ng,0:ng]=raw_barray[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z[ng,0:ng]=Z[0,0:ng]
Z[0:ng+1,ng]=Z[0:ng+1,0]

zmax=np.amax(Z)
dz=contint(0.0,zmax)
jmin=0
jmax=int(zmax/dz)
clevels=np.linspace(0.0,dz*float(jmax),jmax-jmin+1)

im4=ax4.imshow(Z,cmap=cm.seismic,vmin=0.0,vmax=zmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax4)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im4, cax=cax, ticks=clevels)

# -----------------------------------------------------------
t_def=0.4*tmax
t_in = input('  Next time to show (default '+str(t_def)+')? ')
t2 = float(t_in or t_def)
frame=int((t2+0.0001)/dt)

Z=np.empty([ng+1,ng+1])
Z[0:ng,0:ng]=raw_uarray[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z[ng,0:ng]=Z[0,0:ng]
Z[0:ng+1,ng]=Z[0:ng+1,0]

zmax=np.amax(Z)
dz=contint(0.0,zmax)
jmin=0
jmax=int(zmax/dz)
clevels=np.linspace(0.0,dz*float(jmax),jmax-jmin+1)

im2=ax2.imshow(Z,cmap=cm.seismic,vmin=0.0,vmax=zmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im2, cax=cax, ticks=clevels)

Z=np.empty([ng+1,ng+1])
Z[0:ng,0:ng]=raw_barray[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z[ng,0:ng]=Z[0,0:ng]
Z[0:ng+1,ng]=Z[0:ng+1,0]

zmax=np.amax(Z)
dz=contint(0.0,zmax)
jmin=0
jmax=int(zmax/dz)
clevels=np.linspace(0.0,dz*float(jmax),jmax-jmin+1)

im5=ax5.imshow(Z,cmap=cm.seismic,vmin=0.0,vmax=zmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax5)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im5, cax=cax, ticks=clevels)

# -----------------------------------------------------------
t_def=tmax
t_in = input('  Last time to show (default '+str(t_def)+')? ')
t3 = float(t_in or t_def)
frame=int((t3+0.0001)/dt)

Z=np.empty([ng+1,ng+1])
Z[0:ng,0:ng]=raw_uarray[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z[ng,0:ng]=Z[0,0:ng]
Z[0:ng+1,ng]=Z[0:ng+1,0]

zmax=np.amax(Z)
dz=contint(0.0,zmax)
jmin=0
jmax=int(zmax/dz)
clevels=np.linspace(0.0,dz*float(jmax),jmax-jmin+1)

im3=ax3.imshow(Z,cmap=cm.seismic,vmin=0.0,vmax=zmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im3, cax=cax, ticks=clevels)

Z=np.empty([ng+1,ng+1])
Z[0:ng,0:ng]=raw_barray[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z[ng,0:ng]=Z[0,0:ng]
Z[0:ng+1,ng]=Z[0:ng+1,0]

zmax=np.amax(Z)
dz=contint(0.0,zmax)
jmin=0
jmax=int(zmax/dz)
clevels=np.linspace(0.0,dz*float(jmax),jmax-jmin+1)

im6=ax6.imshow(Z,cmap=cm.seismic,vmin=0.0,vmax=zmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax6)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im6, cax=cax, ticks=clevels)

# Add labels for each simulation:
ax1.set_title('$t = {x:.1f}$'.format(x=t1), fontsize=32)
ax2.set_title('$t = {x:.1f}$'.format(x=t2), fontsize=32)
ax3.set_title('$t = {x:.1f}$'.format(x=t3), fontsize=32)
ax1.set_ylabel(r'$\|{\mathbf{u}}\|^2$', fontsize=32)
ax4.set_ylabel(r'$\|{\mathbf{B}}\|^2$', fontsize=32)

ax1.label_outer()
ax2.label_outer()
ax3.label_outer()
ax4.label_outer()
ax5.label_outer()
ax6.label_outer()

# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)
plt.setp(ax6.get_yticklabels(), visible=False)

#fig.subplots_adjust(wspace=0.8, hspace=0.5, top=0.6, bottom=0.05)
fig.subplots_adjust(wspace=0.8, hspace=-0.1)

#=========================================================================
# Save image:
fig.savefig('usqbsq.eps', format='eps', dpi=300)

print()
print(' To view the image, type')
print()
print(' gv usqbsq.eps &')
print()
