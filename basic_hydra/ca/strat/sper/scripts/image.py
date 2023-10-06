#!/usr/bin/env python

# This script plots data in evolution/bb.r4 or zz.r4

#  @@@@   Run from the current job directory   @@@@

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
mpl.rcParams['axes.linewidth'] = 1

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

# Select data to compare:
field_list=['b','zeta']
field_acro=['bb','zz']

print()
print(' Which field do you wish to image')
print()
print(' (1) b')
print(' (2) zeta')
print()
op_in = input(' Option (default 1)? ')
option = int(op_in or 1)

field=field_acro[option-1]
fname='evolution/'+field+'.r4'

print()
t_in = input(' Time to show (default 900)? ')
t = float(t_in or 900.0)

#-----------------------------------------------------------------
# Work out grid resolution (ng) by reading it from parameters.f90:
in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: nx=' in line:
      pline=line

line=pline.split("=")[1]
nx=int(line.split(",")[0])
ny=int(pline.split("=")[2])+1

#-----------------------------------------------------------------
# Open ene.asc file in one directory to get time between frames:
in_file=open('evolution/ecomp.asc','r')
time, ekin, epot, etot = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
nt=len(time)

dt=time[-1]/float(nt-1)
# Frame corresponding to time chosen:
frame=int(t/dt+0.5)
print()
print(' Plotting time',time[frame])

#-----------------------------------------------------------------
# Read data into array for plotting:
in_file=open(fname,'r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()
Z=np.empty([nx,ny])
N=nx*ny
Z=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(nx,ny)

# Work out the overall min/max values:
zmin=np.amin(Z)
zmax=np.amax(Z)
print()
print(' Min & max field levels = '+str(zmin)+', '+str(zmax))
T0=300.0
g=9.80665
dTmax=15.
bmin=-g*dTmax/T0
zm_in = input(' Min level to show (default '+str(bmin)+')? ')
zmin = float(zm_in or bmin)
zm_in = input(' Max level to show (default 0)? ')
zmax = float(zm_in or 0.0)

# Obtain contour levels for plotting the colorbars:
dz=contint(zmin,zmax)
jmin=-int(-zmin/dz)
jmax= int( zmax/dz)
clevels=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

print()
f_in = input(' Starting fraction of x domain to plot (default 0.5)? ')
fx1 = float(f_in or 0.5)
f_in = input(' Ending fraction of x domain to plot (default 0.8515625)? ')
fx2 = float(f_in or 0.8515625)
f_in = input(' Starting fraction of y domain to plot (default 0)? ')
fy1 = float(f_in or 0.0)
f_in = input(' Ending fraction of y domain to plot (default 1)? ')
fy2 = float(f_in or 1.0)

ix1=int(fx1*float(nx)+0.5)
ix2=min(int(fx2*float(nx)+0.5),nx-1)
iy1=int(fy1*float(ny-1)+0.5)
iy2=int(fy2*float(ny-1)+0.5)
print()
print(' Range of x grid points:',ix1,'to',ix2)
print(' Range of y grid points:',iy1,'to',iy2)

#==============================================================================
# Set up figure:
aspect=float(iy2-iy1-1)/float(ix2-ix1-1)
fig1 = plt.figure(1,figsize=[10,10*aspect])
ax1 = fig1.add_subplot(111)
fx1=fx1-0.5
fx2=fx2-0.5
ax1.set_xlim([fx1,fx2])
ax1.set_ylim([fy1,fy2])

ax1.set_xlabel('$x/L_x$', fontsize=20)
ax1.set_ylabel('$y/L_y$', fontsize=20)

# Plot the image in an array with an optional colourbar:
im1=ax1.imshow(Z[ix1:ix2,iy1:iy2].T,cmap=cm.hot,vmin=zmin,vmax=zmax,extent=(fx1,fx2,fy1,fy2),origin='lower',interpolation='bilinear',aspect='auto')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="4%", pad=0.1)
cbar=fig1.colorbar(im1, cax=cax, ticks=clevels)

#=========================================================================
# Save image:
figname=field+'.eps'
fig1.savefig(figname, format='eps', dpi=300)

print()
print(' To view the image, type')
print()
print(' gv '+figname+' &')
print()
