#!/usr/bin/env python

# This script plots data in hovu.r8, containing u_bar(t,phi).

#     @@@@   Run from the current job directory   @@@@

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
# Work out grid resolution (ng) by reading it from parameters.f90:
in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: ng=' in line:
      pline=line

line=pline.split("=")[1]
ng=int(line.split(",")[0])

# Work out number of times in hovmoller plot from size of data:
file_bytes = os.path.getsize('hovu.r8')
nt=int((file_bytes-8)/(8*ng))
print(' Number of times in data = ',nt)

# Read data into array for plotting:
in_file=open('hovu.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z=np.empty([nt,ng])
Z=raw_array[1:nt*ng+1].reshape(nt,ng)

# Work out the overall min/max values:
zmin=np.amin(Z)
zmax=np.amax(Z)

# Centre colourmap at 0:
zmag=max(abs(zmin),zmax)
zmin=-zmag
zmax= zmag

# Obtain contour levels for plotting the colorbars:
dz=contint(zmin,zmax)
jmin=-int(-zmin/dz)
jmax= int( zmax/dz)
clevels=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

# Obtain final time (tmax) in data:
in_file=open('evolution/ecomp.asc','r')
time, dum1, dum2, dum3, dum4=np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

tmax=time[-1]
print(' t_max = ',tmax)

#==============================================================================
# Set up figure:
fig1 = plt.figure(1,figsize=[8,4])
ax1 = fig1.add_subplot(111)
ax1.set_xlim([0.0,tmax])
ax1.set_ylim([-np.pi/2.0,np.pi/2.0])

#ax1.set_xlabel('$t$', fontsize=20)
#ax1.set_ylabel('$\\phi$', fontsize=20)

# Customise tick values:
ax1.yaxis.set_ticks([-np.pi/2.0,-np.pi/4.0,0.0,np.pi/4.0,np.pi/2.0])
ax1.set_yticklabels([r'$-\pi/2$',r'$-\pi/4$',r'$0$',r'$\pi/4$',r'$\pi/2$'],fontsize=20)

# Plot the image in an array with an optional colourbar:
im1=ax1.imshow(Z.T,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(0.0,tmax,-np.pi/2.0,np.pi/2.0),origin='lower',interpolation='bilinear',aspect='auto')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="4%", pad=0.1)
cbar=fig1.colorbar(im1, cax=cax, ticks=clevels)

#=========================================================================
# Save image:
fig1.savefig('hovu.eps', format='eps', dpi=300)

print()
print(' To view the image, type')
print()
print(' gv hovu.eps &')
print()
