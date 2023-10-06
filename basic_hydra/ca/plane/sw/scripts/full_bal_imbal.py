#!/usr/bin/env python

# This script plots the fields for a selected field, and either the
# full field, the balanced field, or the imbalanced field.  Compares
# results in 4 directories, indicated below (2 in SW and 2 in GN),
# at a selected resolution.

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
# Select data to compare:
print()
print(' The following fields may be imaged (full, balanced, imbalanced):')
print()
print(' (1) h;')
print(' (2) zeta;')
print(' (3) delta;')
print(' (4) gamma.')
print()
k=int(input(' Option (default 3)? ') or 3)
k=k-1

field_list=['h','zeta','delta','gamma']
field_acro=['hh','zz','dd','gg']
field=field_list[k]
acron=field_acro[k]

print()
t=float(input(' Time to show (default 500)? ') or 500.0)

print()
ng=int(input(' Resolution (default 256)? ') or 256)

# Define grid:
#xg=np.linspace(-np.pi,np.pi,ng+1)
#yg=xg
#X,Y=np.meshgrid(xg,yg)
N=ng*ng

#=================================================================
# Open ecomp.asc file to get time between frames:
in_file=open('evolution/ecomp.asc','r')
time, e1, e2, e3, e4, e5 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

dt=time[1]-time[0]
# Frame corresponding to time chosen:
frame=int((t+0.0001)/dt)

# Create unique plot output file name:
outfile=field+'_n'+str(ng)+'_t'+str(int(t+0.01))+'.eps'

#=================================================================
# Read data into arrays for plotting:

# Read full field:
in_file=open('evolution/'+acron+'.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z1=np.empty([ng+1,ng+1])
Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

# Read balanced field:
in_file=open('evolution/b'+acron+'.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z2=np.empty([ng+1,ng+1])
Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng)

# Add periodic edges:
Z1[ng,0:ng]=Z1[0,0:ng]
Z1[0:ng+1,ng]=Z1[0:ng+1,0]
Z2[ng,0:ng]=Z2[0,0:ng]
Z2[0:ng+1,ng]=Z2[0:ng+1,0]

# Define imbalanced field:
Z3=Z1-Z2

#====================================================================
# Work out the overall min/max values:
zmin1=np.amin(Z1)
zmin2=np.amin(Z2)
zmin3=np.amin(Z3)
zmax1=np.amax(Z1)
zmax2=np.amax(Z2)
zmax3=np.amax(Z3)

print(zmin1,zmax1)
print(zmin2,zmax2)
print(zmin3,zmax3)

# For full and balanced fields, use same limits:
zmin1=min(zmin1,zmin2)
zmax1=max(zmax1,zmax2)
zmag=max(abs(zmin1),zmax1)
zmin1=-zmag
zmax1= zmag
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

#==============================================================================
# Set up figure:
fig, (ax1, ax2, ax3) = plt.subplots(figsize=[15,6.6], nrows=1, ncols=3)

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

# Plot the images in an array with individual colourbars:
im1=ax1.imshow(Z1,cmap=cm.seismic,vmin=zmin1,vmax=zmax1,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im1, cax=cax, ticks=clevels1)
setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=16)

im2=ax2.imshow(Z2,cmap=cm.seismic,vmin=zmin2,vmax=zmax2,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im2, cax=cax, ticks=clevels2)
setp(cbar.ax.yaxis.set_ticklabels(clevels2), fontsize=16)

im3=ax3.imshow(Z3,cmap=cm.seismic,vmin=zmin3,vmax=zmax3,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im3, cax=cax, ticks=clevels3)
setp(cbar.ax.yaxis.set_ticklabels(clevels3), fontsize=16)

ax1.label_outer()
ax2.label_outer()
ax3.label_outer()

# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
#plt.setp(ax1.get_xticklabels(), visible=False)
#plt.setp(ax2.get_xticklabels(), visible=False)

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
