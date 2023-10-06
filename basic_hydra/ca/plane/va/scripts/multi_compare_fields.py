#!/usr/bin/env python

# This script plots a chosen field (either full, balanced or imbalanced)
# for data at three times in 2 separate directories, as indicated below.

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
# Specify the data directories (need final /):
swdir='/home/dgd/data/hydra/ca/plane/sw/caps/ranpv/ng256kd6dr10r001/'
vadir='/home/dgd/data/hydra/ca/plane/va/caps/ranpv/ng256H0.4kd6dr10r001/'
# Value of H in VA:
hbar=0.4

# Specify grid resolution (ng x ng):
ng=256
N=ng*ng

#=================================================================
# Select data to compare:
print()
print(' Which field do you wish to show:')
print()
print(' (1) h_tilde')
print(' (2) zeta')
print(' (3) delta')
print(' (4) gamma')
print(' (5) gamma_tilde')
print(' (6) q or')
print(' (7) P')

print()
k_in = input(' Option (default 3)? ')
k = int(k_in or 3)

if k==1:
   dfile='hh.r8'
   field='h'
elif k==2:
   dfile='zz.r8'
   field='z'
elif k==3:
   dfile='dd.r8'
   field='d'
elif k==4:
   dfile='gg.r8'
   field='g'
elif k==5:
   dfile='gt.r8'
   field='gt'
elif k==6:
   dfile='qq.r8'
   field='q'
elif k==7:
   dfile='pp.r8'
   field='p'
else:
   print(' *** Not a valid option!  Stopping! ***')
   exit

# Select full, balanced or imbalanced fields:
opt_in = input(' show (1) full, (2) balanced or (3) imbalanced field (default 1)? ')
option = int(opt_in or 1)

if option==2:
   dfile='b'+dfile
   field='b'+field
elif option==3:
   field='i'+field

# Select times to show:
t1_in = input(' 1st time to show (default 50)? ')
t1 = float(t1_in or 50.0)

t2_in = input(' 2nd time to show (default 200)? ')
t2 = float(t2_in or 200.0)

t3_in = input(' 3rd time to show (default 500)? ')
t3 = float(t3_in or 500.0)

# Output filename:
outfile=field+'_n'+str(ng)+'.eps'

#=================================================================
# Open ecomp.asc file in one directory to get time between frames:
in_file=open(swdir+'evolution/ecomp.asc','r')
time, e1, e2, e3, e4, e5 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

dt=time[1]-time[0]
# Frames corresponding to time chosen:
frame1=int((t1+0.0001)/dt)
frame2=int((t2+0.0001)/dt)
frame3=int((t3+0.0001)/dt)

#==============================================================================
# Set up figure:
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(figsize=[16,9.6], nrows=2, ncols=3)

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
# Read data into arrays for plotting:
in_file=open(swdir+'evolution/'+dfile,'r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z1=np.empty([ng+1,ng+1])
Z1[0:ng,0:ng]=raw_array[frame1*(N+1)+1:(frame1+1)*(N+1)].reshape(ng,ng).T
Z2=np.empty([ng+1,ng+1])
Z2[0:ng,0:ng]=raw_array[frame2*(N+1)+1:(frame2+1)*(N+1)].reshape(ng,ng).T
Z3=np.empty([ng+1,ng+1])
Z3[0:ng,0:ng]=raw_array[frame3*(N+1)+1:(frame3+1)*(N+1)].reshape(ng,ng).T

in_file=open(vadir+'evolution/'+dfile,'r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z4=np.empty([ng+1,ng+1])
Z4[0:ng,0:ng]=raw_array[frame1*(N+1)+1:(frame1+1)*(N+1)].reshape(ng,ng).T
Z5=np.empty([ng+1,ng+1])
Z5[0:ng,0:ng]=raw_array[frame2*(N+1)+1:(frame2+1)*(N+1)].reshape(ng,ng).T
Z6=np.empty([ng+1,ng+1])
Z6[0:ng,0:ng]=raw_array[frame3*(N+1)+1:(frame3+1)*(N+1)].reshape(ng,ng).T

if option==3:
   in_file=open(swdir+'evolution/b'+dfile,'r')
   raw_array=np.fromfile(in_file,dtype=np.float64)
   in_file.close()
   Z=np.empty([ng+1,ng+1])
   Z[0:ng,0:ng]=raw_array[frame1*(N+1)+1:(frame1+1)*(N+1)].reshape(ng,ng).T
   Z1=Z1-Z
   Z=np.empty([ng+1,ng+1])
   Z[0:ng,0:ng]=raw_array[frame2*(N+1)+1:(frame2+1)*(N+1)].reshape(ng,ng).T
   Z2=Z2-Z
   Z=np.empty([ng+1,ng+1])
   Z[0:ng,0:ng]=raw_array[frame3*(N+1)+1:(frame3+1)*(N+1)].reshape(ng,ng).T
   Z3=Z3-Z

   in_file=open(vadir+'evolution/b'+dfile,'r')
   raw_array=np.fromfile(in_file,dtype=np.float64)
   in_file.close()
   Z=np.empty([ng+1,ng+1])
   Z[0:ng,0:ng]=raw_array[frame1*(N+1)+1:(frame1+1)*(N+1)].reshape(ng,ng).T
   Z4=Z4-Z
   Z=np.empty([ng+1,ng+1])
   Z[0:ng,0:ng]=raw_array[frame2*(N+1)+1:(frame2+1)*(N+1)].reshape(ng,ng).T
   Z5=Z5-Z
   Z=np.empty([ng+1,ng+1])
   Z[0:ng,0:ng]=raw_array[frame3*(N+1)+1:(frame3+1)*(N+1)].reshape(ng,ng).T
   Z6=Z6-Z

# Add periodic edges:
Z1[ng,0:ng]=Z1[0,0:ng]
Z1[0:ng+1,ng]=Z1[0:ng+1,0]
Z2[ng,0:ng]=Z2[0,0:ng]
Z2[0:ng+1,ng]=Z2[0:ng+1,0]
Z3[ng,0:ng]=Z3[0,0:ng]
Z3[0:ng+1,ng]=Z3[0:ng+1,0]
Z4[ng,0:ng]=Z4[0,0:ng]
Z4[0:ng+1,ng]=Z4[0:ng+1,0]
Z5[ng,0:ng]=Z5[0,0:ng]
Z5[0:ng+1,ng]=Z5[0:ng+1,0]
Z6[ng,0:ng]=Z6[0,0:ng]
Z6[0:ng+1,ng]=Z6[0:ng+1,0]

# Work out the overall min/max values:
zmin1=np.amin(Z1)
zmin2=np.amin(Z2)
zmin3=np.amin(Z3)
zmin4=np.amin(Z4)
zmin5=np.amin(Z5)
zmin6=np.amin(Z6)
zmax1=np.amax(Z1)
zmax2=np.amax(Z2)
zmax3=np.amax(Z3)
zmax4=np.amax(Z4)
zmax5=np.amax(Z5)
zmax6=np.amax(Z6)

# Use different limits for each frame:
zmag=max(abs(zmin1),zmax1)
zmin1=-zmag
zmax1= zmag

zmag=max(abs(zmin2),zmax2)
zmin2=-zmag
zmax2= zmag

zmag=max(abs(zmin3),zmax3)
zmin3=-zmag
zmax3= zmag

zmag=max(abs(zmin4),zmax4)
zmin4=-zmag
zmax4= zmag

zmag=max(abs(zmin5),zmax5)
zmin5=-zmag
zmax5= zmag

zmag=max(abs(zmin6),zmax6)
zmin6=-zmag
zmax6= zmag

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

dz=contint(zmin4,zmax4)
jmin=-int(-zmin4/dz)
jmax= int( zmax4/dz)
clevels4=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin5,zmax5)
jmin=-int(-zmin5/dz)
jmax= int( zmax5/dz)
clevels5=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin6,zmax6)
jmin=-int(-zmin6/dz)
jmax= int( zmax6/dz)
clevels6=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

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

im4=ax4.imshow(Z4,cmap=cm.seismic,vmin=zmin4,vmax=zmax4,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax4)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im4, cax=cax, ticks=clevels4)
#setp(cbar.ax.yaxis.set_ticklabels(clevels4), fontsize=16)

im5=ax5.imshow(Z5,cmap=cm.seismic,vmin=zmin5,vmax=zmax5,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax5)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im5, cax=cax, ticks=clevels5)
#setp(cbar.ax.yaxis.set_ticklabels(clevels5), fontsize=16)

im6=ax6.imshow(Z6,cmap=cm.seismic,vmin=zmin6,vmax=zmax6,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax6)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im6, cax=cax, ticks=clevels6)
#setp(cbar.ax.yaxis.set_ticklabels(clevels6), fontsize=16)

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

# Add labels:
ax1.set_title('$t = {x:.0f}$'.format(x=t1), fontsize=30)
ax2.set_title('$t = {x:.0f}$'.format(x=t2), fontsize=30)
ax3.set_title('$t = {x:.0f}$'.format(x=t3), fontsize=30)
ax1.set_ylabel('$H = 0$', fontsize=30)
ax4.set_ylabel('$H = {x:.1f}$'.format(x=hbar), fontsize=30)

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
