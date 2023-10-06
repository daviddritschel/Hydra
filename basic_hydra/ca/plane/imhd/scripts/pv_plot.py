#!/usr/bin/env python

# This script plots the PV field at a selected time or frame
# from fine-grid data created using genfg.

#========== Perform the generic imports =========
import sys,os,warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.colors as clrs
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

#===========================================================================
print
nx_in = input(' Resolution of ultra-fine grid (default 8192)? ')
nx = int(nx_in or 8192)
dpi_in = input(' DPI to save plot with (default 600)? ')
dpi = int(dpi_in or 600)

frame_in = input(' Frame to show (1, 2, ...; default 1)? ')
frame = int(frame_in or 1)
if frame < 10:
   dataset='qq00'+str(frame)+'.r4'
elif frame < 100:
   dataset='qq0'+str(frame)+'.r4'
else:
   dataset='qq'+str(frame)+'.r4'

opt_in = input(' Add a colourbar (default n)? ')
cbopt=str(opt_in or 'n')

# Read data:
in_file=open('fine/'+dataset,'r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

# Get size to work out y dimension:
file_bytes = os.path.getsize('fine/'+dataset)
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

# Get aspect ratio of the plot:
aspect=yam/np.pi

# Image width (an extra pixel is included for x periodicity):
width=float(8193)/float(600)

# Set up figure:
fig = plt.figure(1,figsize=[width,width])
ax = fig.add_subplot(111)

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

# Plot the image (optionally with a colourbar):
if cbopt == 'n':
   ax.imshow(Z.T,cmap=cm.terrain,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-yam,yam),origin='lower',interpolation='bilinear')
else:
   im=ax.imshow(Z.T,cmap=cm.terrain,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-yam,yam),origin='lower',interpolation='bilinear')
   # Obtain contour levels for plotting the colorbar:
   dz=contint(zmin,zmax)
   jmin=-int(-zmin/dz)
   jmax= int( zmax/dz)
   clevels=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)
   # Add a colourbar with nice intervals:
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("bottom", size="7%", pad=0.1)
   cbar=fig.colorbar(im, cax=cax, ticks=clevels, orientation='horizontal')
   setp(cbar.ax.xaxis.set_ticklabels(clevels), fontsize=32)

#=========================================================================
# Save image:
obase='fine/pv'+str(frame)
fig.savefig(obase+'.eps', format='eps', dpi=dpi, bbox_inches = 'tight')

print()
print(' To view the image, type')
print()
print(' gv '+obase+'.eps &')
print()
print(' It may be worth converting first using')
print()
print(' epstopdf '+obase+'.eps')
print(' ev '+obase+'.pdf &')

