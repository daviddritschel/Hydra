#!/usr/bin/env python

# This script plots <q^2> and <j^2> from data in evolution/norms.asc

#========== Perform the generic imports =========
import warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
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
# Set up figure:
fig = plt.figure(1,figsize=[6.6,6])
ax1 = plt.axes([0.2, 0.2, 0.75, 0.75])

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('$\\langle{q^2}\\rangle,~\\langle{j^2}\\rangle$', fontsize=30)

in_file=open('evolution/norms.asc','r')
time, y1, y2, y3 = np.loadtxt(in_file,dtype=float,unpack=True)
# y3 contains mean-square A; this is much small than mean-square q or j
in_file.close()
xmax_def=time[-1]

xmax_in = input(' Maximum time to show (default '+str(xmax_def)+')? ')
xmax = int(xmax_in or xmax_def)

# Plot:
ax1.plot(time,y1,c='b',lw=2,label='$\\langle{q^2}\\rangle$')
ax1.plot(time,y2,c='r',lw=2,label='$\\langle{j^2}\\rangle$')
ax1.legend(loc='lower left',prop={'size':20})

#=========================================================================
# Determine nice x & y ranges for plotting:
ax1.set_xlim(0.0,xmax)
ax1.set_yscale('log')
ax1.set_box_aspect(1)

# Save figure:
fig.savefig('norms.eps', format='eps', dpi=1200)

print()
print(' To view the image, type')
print()
print(' gv norms.eps &')
print()
