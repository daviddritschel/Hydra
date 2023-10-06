#!/usr/bin/env python

# This script plots the PV field at selected times (specified below
# in the datasets) from images created using genfg.

#========== Perform the generic imports =========
import sys,os,warnings
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

    mpow=0
    rmult=fmax-fmin
    while rmult < 10.0:
       mpow+=1
       rmult=rmult*10.0

    while rmult >= 100.0:
       mpow-=1
       rmult=rmult/10.0

    emag=10.0**(float(-mpow))

    kmult=int(rmult/10.0)

    if kmult < 1:
       ci=emag
    elif kmult < 2:
       ci=2.0*emag
    elif kmult < 4:
       ci=4.0*emag
    elif kmult < 8:
       ci=10.0*emag
    else:
       ci=20.0*emag

    return ci

#=================================================================
print()
nx_in=input(' Resolution of ultra-fine grid (default 8192)? ')
nx=int(nx_in or 8192)

#==============================================================================
# Set up figure:
fig, (ax1, ax2, ax3, ax4) = plt.subplots(figsize=[10,10], nrows=4, ncols=1)

ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.xaxis.set_ticks_position('none')
ax1.yaxis.set_ticks_position('none')

ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.xaxis.set_ticks_position('none')
ax2.yaxis.set_ticks_position('none')

ax3.set_xticklabels([])
ax3.set_yticklabels([])
ax3.xaxis.set_ticks_position('none')
ax3.yaxis.set_ticks_position('none')

ax4.set_xticklabels([])
ax4.set_yticklabels([])
ax4.xaxis.set_ticks_position('none')
ax4.yaxis.set_ticks_position('none')

#=================================================================
# Plot each dataset (top to bottom):

#-----------------
dataset='qq001.r4'

# Read data:
in_file=open(dataset,'r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

# Get size to work out y dimension:
file_bytes = os.path.getsize(dataset)
ny = int((file_bytes-4)/(4*nx))

# Store in 2D array and add periodic edge at x = pi:
Z=np.empty((nx+1,ny))
Z[0:nx,0:ny]=raw_array[1:nx*ny+1].reshape(nx,ny)
Z[nx,0:ny]=Z[0,0:ny]

# Get min/max field values:
zmin=np.amin(Z)
zmax=np.amax(Z)

# Get |y|_max:
yam=np.pi*float(ny)/float(nx)

# Plot the image:
ax1.imshow(Z.T,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-yam,yam),origin='lower',interpolation='bilinear')

#-----------------
dataset='qq003.r4'

# Read data:
in_file=open(dataset,'r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

# Get size to work out y dimension:
file_bytes = os.path.getsize(dataset)
ny = int((file_bytes-4)/(4*nx))

# Store in 2D array and add periodic edge at x = pi:
Z=np.empty((nx+1,ny))
Z[0:nx,0:ny]=raw_array[1:nx*ny+1].reshape(nx,ny)
Z[nx,0:ny]=Z[0,0:ny]

# Get min/max field values:
zmin=np.amin(Z)
zmax=np.amax(Z)

# Get |y|_max:
yam=np.pi*float(ny)/float(nx)

# Plot the image:
ax2.imshow(Z.T,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-yam,yam),origin='lower',interpolation='bilinear')

#-----------------
dataset='qq007.r4'

# Read data:
in_file=open(dataset,'r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

# Get size to work out y dimension:
file_bytes = os.path.getsize(dataset)
ny = int((file_bytes-4)/(4*nx))

# Store in 2D array and add periodic edge at x = pi:
Z=np.empty((nx+1,ny))
Z[0:nx,0:ny]=raw_array[1:nx*ny+1].reshape(nx,ny)
Z[nx,0:ny]=Z[0,0:ny]

# Get min/max field values:
zmin=np.amin(Z)
zmax=np.amax(Z)

# Get |y|_max:
yam=np.pi*float(ny)/float(nx)

# Plot the image:
ax3.imshow(Z.T,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-yam,yam),origin='lower',interpolation='bilinear')

#-----------------
dataset='qq011.r4'

# Read data:
in_file=open(dataset,'r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

# Get size to work out y dimension:
file_bytes = os.path.getsize(dataset)
ny = int((file_bytes-4)/(4*nx))

# Store in 2D array and add periodic edge at x = pi:
Z=np.empty((nx+1,ny))
Z[0:nx,0:ny]=raw_array[1:nx*ny+1].reshape(nx,ny)
Z[nx,0:ny]=Z[0,0:ny]

# Get min/max field values:
zmin=np.amin(Z)
zmax=np.amax(Z)

# Get |y|_max:
yam=np.pi*float(ny)/float(nx)

# Plot the image:
im=ax4.imshow(Z.T,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-yam,yam),origin='lower',interpolation='bilinear')

#-------------------------------------------------
# Obtain contour levels for plotting the colorbar:
dz=contint(zmin,zmax)
jmin=-int(-zmin/dz)
jmax= int( zmax/dz)
clevels=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

# Add a colourbar with nice intervals:
divider = make_axes_locatable(ax4)
cax = divider.append_axes("bottom", size="7%", pad=0.1)
cbar=fig.colorbar(im, cax=cax, ticks=clevels, orientation='horizontal')
setp(cbar.ax.xaxis.set_ticklabels(clevels), fontsize=16)

#=========================================================================
# Save image:
fig.savefig('pv-evo.eps', format='eps', dpi=300)

print()
print( ' To view the image, type')
print()
print( ' gv pv-evo.eps')
print()
