#!/usr/bin/env python

# This script plots the kinetic, potential and total energy versus time,
# optionally comparing with the zonal parts.  Separately plots the
# angular momentum versus time.

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
# Read energies:
in_file=open('evolution/ecomp.asc','r')
t, ekin, epot, etot, angm = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

opt_in = input(' Plot also the zonal counterparts? (default n)? ')
opt = str(opt_in or 'n')

if opt != 'n':
   in_file=open('evolution/zecomp.asc','r')
   tz, ekinz, epotz, etotz, angmz = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

#=================================================================
# Set up figure:
fig, (ax1, ax2) = plt.subplots(figsize=[14.0,6.0], ncols=2)

tmax=float(int(t[-1])+1)

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('$\mathcal{K},~\mathcal{P},~\mathcal{E}$', fontsize=30)
ax1.set_xlim(0.0,tmax)
ax1.plot(t,ekin,c='r',lw=1,label='$\mathcal{K}$')
ax1.plot(t,epot,c='b',lw=1,label='$\mathcal{P}$')
ax1.plot(t,etot,c='k',lw=1,label='$\mathcal{E}$')
ax1.legend(loc='best',prop={'size':20})
ymin,ymax=ax1.get_ylim()
ax1.set_ylim(ymin,ymax)
ax1.set(adjustable='box', aspect=tmax/(ymax-ymin))

ax2.set_xlabel('$t$', fontsize=30)
ax2.set_ylabel('$\mathcal{A}$', fontsize=30)
ax2.set_xlim(0.0,tmax)
ax2.plot(t,angm,c='k',lw=1)
ymin,ymax=ax2.get_ylim()
ax2.set_ylim(ymin,ymax)
ax2.set(adjustable='box', aspect=tmax/(ymax-ymin))

if opt != 'n':
   ax1.plot(tz,ekinz,c='r',lw=1,dashes=[4,4])
   ax1.plot(tz,epotz,c='b',lw=1,dashes=[4,4])
   ax1.plot(tz,etotz,c='k',lw=1,dashes=[4,4])
   ax2.plot(tz,angmz,c='k',lw=1,dashes=[4,4])

# Save image:
fig.savefig('ene_ang.eps', format='eps', dpi=300)

print(' To display the results, type:')
print()
print(' gv ene_ang.eps &')
print()
