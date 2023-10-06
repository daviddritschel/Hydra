#!/usr/bin/env python

# This script plots the divergence field for a selected resolution, 
# from simulations initially balanced or imbalanced.

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
print
print ' Select the initialisation type:'
print
print ' (1) balanced;'
print ' (2) imbalanced.'
print
iopt=int(raw_input(' Option (default 2)? ') or 2)

print
print ' Select the times to show:'
print
if iopt == 1:
   print ' (1) t = 0, 0.5, 1 and 1.5;'
else:
   print ' (1) t = 0.1, 0.5, 1 and 1.5;'
print ' (2) t = 3, 5, 15 and 20.'
print
topt=int(raw_input(' Option (default 1)? ') or 1)

if topt == 1:
   if iopt == 1:
      times=np.array([0.0,0.5,1.0,1.5])
      title=['$t=0$','$t=0.5$','$t=1$','$t=1.5$']
   else:
      times=np.array([0.1,0.5,1.0,1.5])
      title=['$t=0.1$','$t=0.5$','$t=1$','$t=1.5$']
else:
   times=np.array([3.0,5.0,15.0,20.0])
   title=['$t=3$','$t=5$','$t=15$','$t=20$']

print
ng=int(raw_input(' Resolution (default 512)? ') or 512)
print

# Total number of grid points:
N=ng*ng

#=================================================================
# Files containing the data:
if iopt == 1:
   gndata='gn/bal_ng'+str(ng)+'/fine/dd.r8'
   swdata='sw/bal_ng'+str(ng)+'/fine/dd.r8'
else:
   gndata='gn/ng'+str(ng)+'/fine/dd.r8'
   swdata='sw/ng'+str(ng)+'/fine/dd.r8'

# Time between frames:
dt=0.1

# Frames to display:
frames=[]
for t in times:
   frames.append(int((t+1.e-6)/dt))

#=================================================================
# Read in all data to get max abs value for imaging:
d=np.empty((8,ng+1,ng+1))

in_file=open(gndata,'r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
for i in range(4):
   frame=frames[i]
   d[i,0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

in_file=open(swdata,'r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
for i in range(4):
   frame=frames[i]
   d[i+4,0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

# Add periodic edges:
for i in range(8):
   d[i,ng,0:ng]=d[i,0,0:ng]
   d[i,0:ng+1,ng]=d[i,0:ng+1,0]

# Work out the overall max abs value:
dmax = abs(d.max())

# Obtain contour levels for plotting the colorbar:
dd=contint(-dmax,dmax)
jmax=int(dmax/dd)
clevels=np.linspace(-dd*float(jmax),dd*float(jmax),2*jmax+1)

#==============================================================================
# Set up figure:
fig, ax = plt.subplots(figsize=[20,11.3], nrows=2, ncols=4)

ax = ax.flatten()
# Customise tick values and plot images:
for i in range(8):
   ax0=ax[i]
   ax0.set_xticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax0.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=24)
   ax0.set_yticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax0.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=24)
   im=ax0.imshow(d[i],cmap=cm.seismic,vmin=-dmax,vmax=dmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')

# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
for i in range(4):
   ax0=ax[i]
   plt.setp(ax0.get_xticklabels(), visible=False)
   ax0.set_title(title[i], fontsize=30)

# Add information about the code type in the first panels in each row:
if iopt == 1:
   ax0=ax[0]
   ax0.text( 0.1, 0.85, 'GN-bal', transform=ax0.transAxes, fontsize=25, fontname='Times New Roman')
   ax0=ax[4]
   ax0.text( 0.1, 0.85, 'SW-bal', transform=ax0.transAxes, fontsize=25, fontname='Times New Roman')
else:
   ax0=ax[0]
   ax0.text( 0.1, 0.85, 'GN', transform=ax0.transAxes, fontsize=25, fontname='Times New Roman')
   ax0=ax[4]
   ax0.text( 0.1, 0.85, 'SW', transform=ax0.transAxes, fontsize=25, fontname='Times New Roman')

ax0=ax[5]
ax0.xaxis.labelpad = 80
ax0.set_xlabel('.', fontsize=1)
ax0=ax[6]
ax0.xaxis.labelpad = 80
ax0.set_xlabel('.', fontsize=1)

for i in range(3):
   plt.setp(ax[i+1].get_yticklabels(), visible=False)
   plt.setp(ax[i+5].get_yticklabels(), visible=False)

#fig.tight_layout()

fig.subplots_adjust(bottom=0.2, top=1.0, left=0.1, right=0.9,
                    wspace=-0.1, hspace=-0.1)

# Add colourbar underneath:
ax_cbar = fig.add_axes([0.25, 0.05, 0.5, 0.03])
cbar=fig.colorbar(im, cax=ax_cbar, ticks=clevels, orientation='horizontal')
setp(cbar.ax.xaxis.set_ticklabels(clevels), fontsize=24)

#fig.subplots_adjust(wspace=0.8, hspace=0.5, top=0.6, bottom=0.05)
#fig.subplots_adjust(wspace=0.8, hspace=0.1)

#=========================================================================
# Save image:

print
print ' To view the image, type'
print
if topt == 1:
   if iopt == 1:
      fig.savefig('d_bal'+str(ng)+'early.eps', format='eps', dpi=300)
      print ' gv d_bal'+str(ng)+'early.eps'
   else:
      fig.savefig('d'+str(ng)+'early.eps', format='eps', dpi=300)
      print ' gv d'+str(ng)+'early.eps'
else:
   if iopt == 1:
      fig.savefig('d_bal'+str(ng)+'late.eps', format='eps', dpi=300)
      print ' gv d_bal'+str(ng)+'late.eps'
   else:
      fig.savefig('d'+str(ng)+'late.eps', format='eps', dpi=300)
      print ' gv d'+str(ng)+'late.eps'
print
