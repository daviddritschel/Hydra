#!/usr/bin/env python

# This script plots data in ro-fr-hm.asc, with time filtering

#========== Perform the generic imports =========
import numpy as np
from scipy.ndimage.filters import uniform_filter1d
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
# Read diagnostics:
in_file=open('evolution/ro-fr-hm.asc','r')
t, ro, fr, hmin, hmax, zrms, drms = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

tsamp_in = input(' Window size in time (default 10)? ')
tsamp = float(tsamp_in or 10.0)
nsamp = int(tsamp/(t[1]-t[0])+0.5)

# Perform running average:
aro = uniform_filter1d(ro, size=nsamp, mode='nearest')
afr = uniform_filter1d(fr, size=nsamp, mode='nearest')
ahmin = uniform_filter1d(hmin, size=nsamp, mode='nearest')
ahmax = uniform_filter1d(hmax, size=nsamp, mode='nearest')
azrms = uniform_filter1d(zrms, size=nsamp, mode='nearest')
adrms = uniform_filter1d(drms, size=nsamp, mode='nearest')

#=================================================================
# Set up figure:
fig, (ax1, ax2, ax3) = plt.subplots(figsize=[20.0,6.0], ncols=3)

tmax=float(int(t[-1])+1)

ymax=max(np.amax(aro),np.amax(afr))
ymax=0.05*(int(20.0*ymax)+1.0)
ax1.set(adjustable='box', aspect=tmax/ymax)
ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('$\mathrm{Ro},~\mathrm{Fr}$', fontsize=30)
ax1.set_xlim(0.0,tmax)
ax1.set_ylim(0.0,ymax)
ax1.plot(t,aro,c='k',lw=1,label='$\mathrm{Ro}$')
ax1.plot(t,afr,c='b',lw=1,label='$\mathrm{Fr}$')
ax1.legend(loc='lower right',prop={'size':20})

ymax=max(np.amax(abs(ahmin)),np.amax(ahmax))
ymax=0.05*(int(20.0*ymax)+1.0)
ax2.set(adjustable='box', aspect=tmax/(2.0*ymax))
ax2.set_xlabel('$t$', fontsize=30)
ax2.set_ylabel('$h_{\mathsf{min}},~h_{\mathsf{max}}$', fontsize=30)
ax2.set_xlim(0.0,tmax)
ax2.set_ylim(-ymax,ymax)
ax2.plot(t,ahmin,c='k',lw=1)
ax2.plot(t,ahmax,c='k',lw=1)

ymax=max(np.amax(azrms),np.amax(adrms))
ymax=0.1*(int(10.0*ymax)+1.0)
ax3.set(adjustable='box', aspect=tmax/ymax)
ax3.set_xlabel('$t$', fontsize=30)
ax3.set_ylabel('$\\zeta_{\mathsf{rms}},~\\delta_{\mathsf{rms}}$', fontsize=30)
ax3.set_xlim(0.0,tmax)
ax3.set_ylim(0.0,ymax)
ax3.plot(t,azrms,c='k',lw=1,label='$\\zeta_{\mathsf{rms}}$')
ax3.plot(t,adrms,c='b',lw=1,label='$\\delta_{\mathsf{rms}}$')
ax3.legend(loc='lower right',prop={'size':20})

# Save image:
#plt.tight_layout()
#fig.subplots_adjust(wspace=0.1, hspace=0.1)

fig.savefig('diag.eps', format='eps', dpi=300)

print(' To display the results, type:')
print()
print(' gv diag.eps &')
print()
