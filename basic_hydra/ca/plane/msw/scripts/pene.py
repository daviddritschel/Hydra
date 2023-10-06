#!/usr/bin/env python

# This script plots the kinetic, potential, magnetic and total energy
# versus time, with the possibility of comparing data in different
# directories as specified below.

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

# Set up figure:
fig1 = plt.figure(1,figsize=[8,8])
ax1 = fig1.add_subplot(111)

#=================================================================
# For working out maximum energy and time:
emax=0.0
tmax=0.0

# Read data and plot:
for m,dfile in enumerate(dir_list):
   in_file=open(dir_list[m]+'evolution/ecomp.asc','r')
   t, ekin, epot, emag, etot = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   # Scale time:
   t=rossby[m]*t
   tmax=max(tmax,t[-1])

   emax=max(emax,np.amax(etot))

   if m==0:
      ax1.plot(t,ekin,color=colour_list[m],lw=2,linestyle='solid',label='$\mathcal{E}_u$')
      ax1.plot(t,emag,color=colour_list[m],lw=2,linestyle='dotted',label='$\mathcal{E}_b$')
      ax1.plot(t,epot,color=colour_list[m],lw=2,linestyle='dashed',label='$\mathcal{E}_h$')
      ax1.plot(t,etot,color=colour_list[m],lw=2,linestyle='dashdot',label='$\mathcal{E}$')
   else:
      ax1.plot(t,ekin,color=colour_list[m],lw=2,linestyle='solid')
      ax1.plot(t,emag,color=colour_list[m],lw=2,linestyle='dotted')
      ax1.plot(t,epot,color=colour_list[m],lw=2,linestyle='dashed')
      ax1.plot(t,etot,color=colour_list[m],lw=2,linestyle='dashdot')

ax1.set_xlim(0.0,tmax)
emax=1.02*emax
emin=1.e-6*emax
ax1.set_ylim(emin,emax)
ax1.set_yscale('log')

ax1.legend(loc='best',prop={'size':25})
ax1.set(adjustable='box', aspect=tmax/emax)

ax1.set_xlabel('$\\tilde{t}$', fontsize=30)
ax1.set_ylabel('$\mathcal{E}_u,\,\,\,\mathcal{E}_b,\,\,\,\mathcal{E}_h,\,\,\,\mathcal{E}$', fontsize=30)

# Save image:
fig1.savefig('ene.eps', format='eps', dpi=300)

print(' To display the results, type:')
print()
print(' gv ene.eps &')
print()
