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
mpl.rcParams['text.latex.preamble'] = r'\usepackage{bm}'

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
data_files=['evolution/','../ng512kd3eps0.1gamma4.0bzrat0.0/evolution/']
# Corresponding line styles:
dashlist=[(1,0.0001),(8,3),(2,4,2),(1,3,3,1)]
# Corresponding labels for the plots:
label1='$\\nabla\\cdot( h{{\\bm{B}} ) =0}$'
label2='$\\nabla\\cdot( h{{\\bm{B}} ) \\neq{0}}$'

# Set up figure:
fig1 = plt.figure(1,figsize=[8,8])
ax1 = fig1.add_subplot(111)

# For working out maximum energy and time:
emax=0.0
emin=1.e20
tmax=0.0

# Read data and plot:
for m,dfile in enumerate(data_files):
   in_file=open(data_files[m]+'ecomp.asc','r')
   t, ekin, epot, emag, etot = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   emin=min(emin,np.amin(etot))
   emax=max(emax,np.amax(etot))
   tmax=max(tmax,t[-1])

   if m==0:
      ax1.plot(t,etot,c='k',lw=2,dashes=dashlist[m],label=label1)
   else:
      ax1.plot(t,etot,c='k',lw=2,dashes=dashlist[m],label=label2)

ax1.set_xlim(0.0,tmax)
emin=0.98*emin
emax=1.02*emax
ax1.set_ylim(emin,emax)

ax1.legend(loc='best',prop={'size':25})
ax1.set(adjustable='box', aspect=tmax/(emax-emin))

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('$\mathcal{E}$', fontsize=30)

# Save image:
fig1.savefig('tene.eps', format='eps', dpi=300)

print(' To display the results, type:')
print()
print(' gv tene.eps &')
print()
