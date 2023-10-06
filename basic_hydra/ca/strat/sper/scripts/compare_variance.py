#!/usr/bin/env python

# This script plots the variances <b^2> and <zeta^2> comparing 
# selected directories specified below.

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
# Specify the data directories:
dir_list=['1x','2x','4x','8x']
# Corresponding labels:
label_list=['1x','2x','4x','8x']
# Corresponding colours (allow up to 4 directories):
colorlist=['k','b','r','m']

#=================================================================
# Open variance.asc file in one directory to get the final time:
in_file=open(dir_list[0]+'/evolution/variance.asc','r')
time, bbl2, zzl2=np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

tmax=time[-1]

# Set up figure:
fig, (ax1, ax2) = plt.subplots(figsize=[12.5,6], nrows=1, ncols=2)

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('$\\langle{b^2}\\rangle$', fontsize=30)
ax1.set_xlim(0.0,tmax)
#ax1.yaxis.set(major_locator=MultipleLocator(0.1), \
#                 major_formatter=FormatStrFormatter('%1.1f'))
ax1.set_box_aspect(1)

ax2.set_xlabel('$t$', fontsize=30)
ax2.set_ylabel('$\\langle{\\zeta^2}\\rangle$', fontsize=30)
ax2.set_xlim(0.0,tmax)
#ax2.yaxis.set(major_locator=MultipleLocator(0.1), \
#                 major_formatter=FormatStrFormatter('%1.1f'))
ax2.set_box_aspect(1)

#=================================================================
# Loop over directories and plot results:
for m,dir in enumerate(dir_list):
   # Open input file and read data:
   in_file=open(dir+'/evolution/variance.asc','r')
   time, bbl2, zzl2 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   ax1.plot(time,bbl2,c=colorlist[m],lw=2,label=label_list[m])
   ax1.legend(loc='lower left',prop={'size':20})
   ax2.plot(time,zzl2,c=colorlist[m],lw=2)

#=========================================================================
# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0)

fig.savefig('l2.eps',  format='eps', dpi=1200)

print()
print(' To view the image, type')
print()
print(' gv l2.eps &')
print()
