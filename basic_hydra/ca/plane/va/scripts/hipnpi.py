#!/usr/bin/env python

# This script compared h_i, p_n and P_i at a chosen time.

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

#=================================================================
# Specify grid resolution (ng x ng):
print()
ng_in = input(' Resolution (default 256)? ')
ng = int(ng_in or 256)

N=ng*ng

# Select time to show:
t_in = input(' Time to show (default 500)? ')
t = float(t_in or 500.0)

# Output filename:
outfile='hipnpi_n'+str(ng)+'_t'+str(int(t+0.01))+'.eps'

#=================================================================
# Open ecomp.asc file to get time between frames:
in_file=open('evolution/ecomp.asc','r')
time, e1, e2, e3, e4, e5 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

dt=time[1]-time[0]
# Frame corresponding to time chosen:
frame=int((t+0.0001)/dt)

#==============================================================================
# Set up figure:
fig, (ax1, ax2, ax3) = plt.subplots(figsize=[15,5], nrows=1, ncols=3)

# Customise tick values:
ax1.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

ax1.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

#=================================================================
# Read data into arrays for plotting:
Z1=np.empty([ng+1,ng+1])
in_file=open('evolution/hh.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z=np.empty([ng+1,ng+1])
in_file=open('evolution/bhh.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z1=Z1-Z

Z2=np.empty([ng+1,ng+1])
in_file=open('evolution/pn.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

Z3=np.empty([ng+1,ng+1])
in_file=open('evolution/pp.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z=np.empty([ng+1,ng+1])
in_file=open('evolution/bpp.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z3=Z3-Z

# Add periodic edges:
Z1[ng,0:ng]=Z1[0,0:ng]
Z1[0:ng+1,ng]=Z1[0:ng+1,0]
Z2[ng,0:ng]=Z2[0,0:ng]
Z2[0:ng+1,ng]=Z2[0:ng+1,0]
Z3[ng,0:ng]=Z3[0,0:ng]
Z3[0:ng+1,ng]=Z3[0:ng+1,0]

# Work out the overall min/max values:
zmin1=np.amin(Z1)
zmin2=np.amin(Z2)
zmin3=np.amin(Z3)
zmax1=np.amax(Z1)
zmax2=np.amax(Z2)
zmax3=np.amax(Z3)

# Use different limits for each simulation:
zmag=max(abs(zmin1),zmax1)
zmin1=-zmag
zmax1= zmag

zmag=max(abs(zmin2),zmax2)
zmin2=-zmag
zmax2= zmag

zmag=max(abs(zmin3),zmax3)
zmin3=-zmag
zmax3= zmag

# Obtain contour levels for plotting the colorbars:
dz=contint(zmin1,zmax1)
jmin=-int(-zmin1/dz)
jmax= int( zmax1/dz)
clevels1=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin2,zmax2)
jmin=-int(-zmin2/dz)
jmax= int( zmax2/dz)
clevels2=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin3,zmax3)
jmin=-int(-zmin3/dz)
jmax= int( zmax3/dz)
clevels3=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

# Plot the images in an array with individual colourbars:
im1=ax1.imshow(Z1,cmap=cm.seismic,vmin=zmin1,vmax=zmax1,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im1, cax=cax, ticks=clevels1)
#setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=16)

im2=ax2.imshow(Z2,cmap=cm.seismic,vmin=zmin2,vmax=zmax2,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im2, cax=cax, ticks=clevels2)
#setp(cbar.ax.yaxis.set_ticklabels(clevels2), fontsize=16)

im3=ax3.imshow(Z3,cmap=cm.seismic,vmin=zmin3,vmax=zmax3,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im3, cax=cax, ticks=clevels3)
#setp(cbar.ax.yaxis.set_ticklabels(clevels3), fontsize=16)

# Add titles:

# 
# 
# 

ax1.set_title('$\\tilde{h}_{\mathsf{i}}$', fontsize=20)
ax2.set_title('$p_{n}$', fontsize=20)
ax3.set_title('$P_{\mathsf{i}}$', fontsize=20)

ax1.label_outer()
ax2.label_outer()
ax3.label_outer()

# Fine-tune figure; hide y ticks for right plots
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

#fig.subplots_adjust(wspace=0.8, hspace=0.5, top=0.6, bottom=0.05)
fig.subplots_adjust(wspace=0.8, hspace=-0.1)

#=========================================================================
# Save image:
fig.savefig(outfile, format='eps', dpi=300)

print()
print(' To view the image, type')
print()
print(' gv ',outfile)
print()
