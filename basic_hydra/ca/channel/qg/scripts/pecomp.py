#!/usr/bin/env python

# This script plots the kinetic, potential and total energy
# versus time, with the possibility of comparing data in different
# directories as specified below.

#========== Perform the generic imports =========
import warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

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
# Specify data directories here:
dir_list=['nx1024ny256kd1.0beta0.4r001/']
ndir=len(dir_list)

# Corresponding line styles:
dashlist=[(1,0.0001),(8,3),(1,3,3,1)]

# Corresponding labels for the plots:
lab_list=['${\\beta/c}=4$']

# Set up figure:
fig, (ax1, ax2, ax3) = plt.subplots(figsize=[18.2,6], nrows=1, ncols=3)

# For working out maximum energy and time:
emax=0.0
emin=1.e20
epmax=0.0
epmin=1.e20
ekmax=0.0
ekmin=1.e20
tmax=0.0

# Read data and plot:
for m,direc in enumerate(dir_list):
   in_file=open(direc+'evolution/ene.asc','r')
   t, ekin, epot, etot = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   emin=min(emin,np.amin(etot))
   emax=max(emax,np.amax(etot))
   epmin=max(epmin,np.amin(epot))
   epmax=max(epmax,np.amax(epot))
   ekmin=min(ekmin,np.amin(ekin))
   ekmax=max(ekmax,np.amax(ekin))
   tmax=max(tmax,t[-1])

   ax1.plot(t,etot,c='k',lw=2,dashes=dashlist[m])
   ax2.plot(t,epot,c='k',lw=2,dashes=dashlist[m],label=lab_list[m])
   ax3.plot(t,ekin,c='k',lw=2,dashes=dashlist[m])

ax1.set_xlim(0.0,tmax)
ax2.set_xlim(0.0,tmax)
ax3.set_xlim(0.0,tmax)
sfac=0.02
de=sfac*(emax-emin)
emin=emin-de
emax=emax+de
ax1.set_ylim(emin,emax)
de=sfac*(epmax-epmin)
epmax=epmax+de
ax2.set_ylim(epmin,epmax)
de=sfac*(ekmax-ekmin)
ekmin=ekmin-de
ekmax=ekmax+de
ax3.set_ylim(ekmin,ekmax)

ax1.set(adjustable='box', aspect=tmax/(emax-emin))
ax2.legend(loc='best',prop={'size':20})
ax2.set(adjustable='box', aspect=tmax/(epmax-epmin))
ax3.set(adjustable='box', aspect=tmax/(ekmax-ekmin))

ax1.set_xlabel('$\\varepsilon t$', fontsize=30)
ax1.set_ylabel('$\mathcal{E}$', fontsize=30)
ax2.set_xlabel('$\\varepsilon t$', fontsize=30)
ax2.set_ylabel('$\mathcal{P}$', fontsize=30)
ax3.set_xlabel('$\\varepsilon t$', fontsize=30)
ax3.set_ylabel('$\mathcal{K}$', fontsize=30)

# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0)

#fig.savefig('ecomp.eps',  format='eps', dpi=1200)
fig.savefig('ecomp.png')

print(' To display the results, type:')
print()
#print(' gv ecomp.eps &')
print(' eom ecomp.png &')
print()
