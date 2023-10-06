#!/usr/bin/env python

# This script plots the Rossby and Froude numbers, and the min/max values
# of h, comparing data in different directories as specified below.

#========== Perform the generic imports =========
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

# set tick label size:
label_size = 24
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1
# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 1
# set axes width:
mpl.rcParams['axes.linewidth'] = 3

#=================================================================
# Specify the data directories (need final /):
dir_list=['ng512kd2eps0.400gamma2.0bzrat0.0/','ng512kd4eps0.100gamma2.0bzrat0.0/','ng512kd8eps0.025gamma2.0bzrat0.0/']
ndir=len(dir_list)

# Corresponding colours:
colour_list=['k','b','m']

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
# Set up figures:
fig1 = plt.figure(1,figsize=[8,8])
ax1 = fig1.add_subplot(111)

fig2 = plt.figure(2,figsize=[8,8])
ax2 = fig2.add_subplot(111)

# For working out min/max values of the diagnostics:
tmax=0.0
romax=0.0
frmax=0.0
hhmax=0.0
hhmin=0.0

# Read data and plot:
for m,dfile in enumerate(dir_list):
   in_file=open(dir_list[m]+'evolution/ro-fr-hm.asc','r')
   t, ro, fr, hn, hp = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   # Scale time:
   t=rossby[m]*t
   tmax=max(tmax,t[-1])

   romax=max(romax,np.amax(ro))
   frmax=max(frmax,np.amax(fr))
   hhmax=max(hhmax,np.amax(hp))
   hhmin=min(hhmin,np.amin(hn))

   ax1.plot(t,ro,lw=2,color=colour_list[m],label=lab_list[m])
   ax1.plot(t,fr,lw=2,color=colour_list[m],dashes=(8,3))
   ax2.plot(t,hp,lw=2,color=colour_list[m],label=lab_list[m])
   ax2.plot(t,hn,lw=2,color=colour_list[m],dashes=(8,3))

ax1.set_xlim(0.0,tmax)
ymax1=1.02*max(romax,frmax)
ax1.set_ylim(0.0,ymax1)

ax2.set_xlim(0.0,tmax)
ymax2=1.02*hhmax
ymin2=1.02*hhmin
ax2.set_ylim(ymin2,ymax2)

ax1.legend(loc='best',prop={'size':25})
ax1.set(adjustable='box', aspect=tmax/ymax1)

ax2.legend(loc='best',prop={'size':25})
ax2.set(adjustable='box', aspect=tmax/(ymax2-ymin2))

ax1.set_xlabel('$\\tilde{t}$', fontsize=30)
ax1.set_ylabel('$\\rm{Ro},\,\,\,\\rm{Fr}$', fontsize=30)

ax2.set_xlabel('$\\tilde{t}$', fontsize=30)
ax2.set_ylabel('$h_{max},\,\,\,h_{min}$', fontsize=30)

# Save images:
fig1.savefig('rofr.eps', format='eps', dpi=300)
fig2.savefig('hext.eps', format='eps', dpi=300)

print(' To display the results, type:')
print()
print(' gv rofr.eps &')
print(' gv hext.eps &')
print()
