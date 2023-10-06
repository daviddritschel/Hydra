#!/usr/bin/env python

# This script plots a chosen field for data at three times in
# several separate directories --- specified below.

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
mpl.rcParams['text.latex.preamble'] = r'\usepackage{bm}'

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
dir_list=['ng512kd2eps0.400gamma2.0bzrat0.0/','ng512kd4eps0.100gamma2.0bzrat0.0/','ng512kd8eps0.025gamma2.0bzrat0.0/']
ndir=len(dir_list)

# Corresponding Rossby numbers for plotting dimensionless time:
rossby=np.array([0.4,0.1,0.025])

# Corresponding labels for the plots:
lab_list=['$L_D=1/2$','$L_D=1/4$','$L_D=1/8$']

# Specify grid resolution (ng x ng):
print()
ng_in = input(' Resolution (default 512)? ')
ng = int(ng_in or 512)
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
print(' (5) q_l')
print(' (6) j_z')
print(' (7) B_perp = -div(h(B_x,B_y)), or')
print(' (8) B_z at z = 0')

print()
k_in = input(' Option (default 2)? ')
k = int(k_in or 2)

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
   dfile='qq.r8'
   field='ql'
elif k==6:
   dfile='jz.r8'
   field='jz'
elif k==7:
   dfile='bp.r8'
   field='bperp'
elif k==8:
   dfile='bz.r8'
   field='bz'
else:
   print(' *** Not a valid option!  Stopping! ***')
   exit

# Select times to show:
print(' Specify three dimensionless times, Ro*t, to show.')
print()
t1_in = input(' 1st time to show (default 5)? ')
t1 = float(t1_in or 5.0)

t2_in = input(' 2nd time to show (default 10)? ')
t2 = float(t2_in or 10.0)

t3_in = input(' 3rd time to show (default 25)? ')
t3 = float(t3_in or 25.0)
t=np.array([t1,t2,t3])

# Output filename:
outfile=field+'_n'+str(ng)+'.eps'

#=================================================================

# Read data into arrays for plotting:
d=np.empty((3*ndir,ng+1,ng+1))
i=0
for k,direc in enumerate(dir_list):
   # Open ecomp.asc file to get time between frames:
   in_file=open(direc+'evolution/ecomp.asc','r')
   time, e1, e2, e3, e4 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   dt=rossby[k]*(time[1]-time[0])

   # Open main data file and read:
   in_file=open(direc+'evolution/'+dfile,'r')
   raw_array=np.fromfile(in_file,dtype=np.float64)
   in_file.close()

   # Read selected frames:
   for k in range(3):
      j=int((t[k]+0.0001)/dt)
      d[i,0:ng,0:ng]=raw_array[j*(N+1)+1:(j+1)*(N+1)].reshape(ng,ng).T
      if field == 'bz':
         # Remove average:
         davg=np.mean(d[i,:,:])
         print(davg)
         d[i,0:ng,0:ng]=d[i,0:ng,0:ng]-davg
      # Add periodic edges:
      d[i,ng,0:ng]=d[i,0,0:ng]
      d[i,0:ng+1,ng]=d[i,0:ng+1,0]
      i+=1

#==============================================================================
# Set up figure:
fig, ax = plt.subplots(figsize=[16,4.8*ndir], nrows=ndir, ncols=3)
ax = ax.flatten()

# Customise tick values and plot images:
for i in range(3*ndir):
   ax0=ax[i]
   ax0.set_xticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax0.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=24)
   ax0.set_yticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax0.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=24)

   dmax = np.amax(abs(d[i]))
   # Obtain contour levels for plotting the colorbar:
   dd=contint(-dmax,dmax)
   jmax=int(dmax/dd)
   clevels=np.linspace(-dd*float(jmax),dd*float(jmax),2*jmax+1)
   # Plot image:
   im=ax0.imshow(d[i],cmap=cm.seismic,vmin=-dmax,vmax=dmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   # Add individual colorbars:
   divider = make_axes_locatable(ax0)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im, cax=cax, ticks=clevels)
   # Suppress ticks on outer edges:
   ax0.label_outer()

# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
for i in range(3*(ndir-1)):
   ax0=ax[i]
   plt.setp(ax0.get_xticklabels(), visible=False)

# Add times in the titles:
for i in range(3):
   ax0=ax[i]
   ax0.set_title('${{\\tilde{{t}}}} = {x:.0f}$'.format(x=t[i]), fontsize=30)

# Add y labels:
for k in range(ndir):
   ax0=ax[3*k]
   ax0.set_ylabel(lab_list[k], fontsize=30)

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
