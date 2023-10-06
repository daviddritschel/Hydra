#!/usr/bin/env python

# This script plots |B|_max/U_0 and zeta_max/f, using data in
# evolution/norms.asc and in evolution/ro-fr-hm.asc in different
# directories specified below.

# *** Here, U_0 is provided below - change if necessary ***

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
# Specify U_0:
U_0=0.308425137533

# Specify data directories here:
dir_list=['ng1024kd4eps0.100gamma2.0bzrat0.0/','ng512kd4eps0.100gamma2.0bzrat0.0/','ng256kd4eps0.100gamma2.0bzrat0.0/']
ndir=len(dir_list)

# Corresponding line styles:
dashlist=[(1,0.0001),(8,3),(1,3,3,1)]

# Corresponding Rossby numbers for plotting dimensionless time:
rossby=np.array([0.1,0.1,0.1])

# Corresponding labels for the plots:
lab_list=['${\\mathrm Rm}=3200$','${\\mathrm Rm}=800$','${\\mathrm Rm}=200$']

# Set up figure:
fig, (ax1, ax2) = plt.subplots(figsize=[12.3,6], nrows=1, ncols=2)

# For working out min/max values and time:
zmax=0.0
zmin=1.e20
bmax=0.0
bmin=1.e20
tmax=0.0

# Read data and plot:
for m,direc in enumerate(dir_list):

   in_file=open(direc+'evolution/ro-fr-hm.asc','r')
   t, z, d1, d2, d3 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   in_file=open(direc+'evolution/norms.asc','r')
   t, b, d1, d2, d3 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   t=rossby[m]*t
   b=b/U_0
   
   zmin=min(zmin,np.amin(z))
   zmax=max(zmax,np.amax(z))
   bmin=min(bmin,np.amin(b))
   bmax=max(bmax,np.amax(b))
   tmax=max(tmax,t[-1])

   ax1.plot(t,z,c='k',lw=2,dashes=dashlist[m])
   ax2.plot(t,b,c='k',lw=2,dashes=dashlist[m],label=lab_list[m])

ax1.set_xlim(0.0,tmax)
ax2.set_xlim(0.0,tmax)
sfac=0.02
dz=sfac*(zmax-zmin)
zmin=zmin-dz
zmax=zmax+dz
db=sfac*(bmax-bmin)
bmin=bmin-db
bmax=bmax+db
ax1.set_ylim(zmin,zmax)
ax2.set_ylim(bmin,bmax)

ax1.set_yscale('log')

ax1.set(adjustable='box', aspect=tmax/np.log10(zmax/zmin))
#ax1.set(adjustable='box', aspect=tmax/(zmax-zmin))
ax2.set(adjustable='box', aspect=tmax/(bmax-bmin))
ax2.legend(loc='best',prop={'size':20})

ax1.set_xlabel('$\\varepsilon t$', fontsize=30)
ax1.set_ylabel('$\mathcal{\\zeta}_{\\mathsf{max}}/f$', fontsize=30)
ax2.set_xlabel('$\\varepsilon t$', fontsize=30)
ax2.set_ylabel('$\\|{\\bm{B}}\\|_{\\mathsf{max}}/U_0$', fontsize=30)

# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0)

fig.savefig('norms.eps',  format='eps', dpi=1200)

print(' To display the results, type:')
print()
print(' gv norms.eps &')
print()
