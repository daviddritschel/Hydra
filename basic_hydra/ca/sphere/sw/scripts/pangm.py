#======================================================================
# Plots A vs t from ecomp.asc
#======================================================================

#=====perform various generic imports=====
import warnings
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

## global settings

# set tick label size:
label_size = 25
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

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

#===================================================================
#Read data:
dfile = open('ecomp.asc','r')
time, d1, d2, d3, angm = np.loadtxt(dfile,dtype=float,unpack=True)
dfile.close()

#-------------------------------------------------------------------
# Plot results:
fig1 = plt.figure(1,figsize=[12,8])
ax1 = plt.axes([0.2, 0.2, 0.7, 0.7])
#plt.gca().yaxis.set_major_locator( MaxNLocator(nbins = 4) )

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('${\mathcal{A}}$', fontsize=30)
#ax1.set_xlim(t[0],t[-1])
#ax1.set_ylim(0.0,1.04)

ax1.plot(time,angm,c='k',lw=2)

plt.savefig('angm.eps', format='eps', dpi=600, bbox_inches = 'tight')

plt.show()
