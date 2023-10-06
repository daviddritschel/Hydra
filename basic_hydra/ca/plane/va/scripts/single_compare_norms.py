#!/usr/bin/env python

# This script plots rms field norms for a chosen field from data
# in the directories specified below.

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
# Specify the data directories (need final /):
basedir='/home/dgd/data/hydra/ca/plane/'
dir_list=[basedir+'sw/caps/ranpv/ng256kd6dr10r001/', \
          basedir+'va/caps/ranpv/ng256H0.1kd6dr10r001/', \
          basedir+'va/caps/ranpv/ng256H0.2kd6dr10r001/', \
          basedir+'va/caps/ranpv/ng256H0.4kd6dr10r001/']
# Corresponding titles (supplemented by chosen field):
label_list=['$H=0$','$H=0.1$','$H=0.2$','$H=0.4$']
# Corresponding colours (allow up to 4 directories):
colorlist=['k','b','r','m']

#=================================================================
# Select data to plot:
print()
print(' This script compares field norms in different directories.')
print()
print(' Which field do you wish to examine:')
print()
print(' (1) h_tilde')
print(' (2) zeta')
print(' (3) delta')
print(' (4) gamma')
print(' (5) gamma_tilde, or')
print(' (6) P')

print()
opt_in = input(' Option (default 3)? ')
field = int(opt_in or 3)
field = field-1

# Select full, balanced or imbalanced results:
print()
opt_in = input(' show (1) full, (2) balanced or (3) imbalanced norm (default 1)? ')
option = int(opt_in or 1)

opt_in = input(' Add a legend (1 = yes, 0 = no; default 0)? ')
leg = int(opt_in or 0)
   
file_prefix=['h','z','d','g','t','p']

if option==1:
   outfile=file_prefix[field]+'_rms.eps'
   ylabel=['$h_{\mathsf{rms}}$','$\\zeta_{\mathsf{rms}}$','$\\delta_{\mathsf{rms}}$','$\\gamma_{\mathsf{rms}}$','$\\tilde\\gamma_{\mathsf{rms}}$','$P_{\mathsf{rms}}$']
elif option==2:
   outfile='b'+file_prefix[field]+'_rms.eps'
   ylabel=['$h_{\mathsf{b,rms}}$','$\\zeta_{\mathsf{b,rms}}$','$\\delta_{\mathsf{b,rms}}$','$\\gamma_{\mathsf{b,rms}}$','$\\tilde\\gamma_{\mathsf{b,rms}}$','$P_{\mathsf{b,rms}}$']
else:
   outfile='i'+file_prefix[field]+'_rms.eps'
   ylabel=['$h_{\mathsf{i,rms}}$','$\\zeta_{\mathsf{i,rms}}$','$\\delta_{\mathsf{i,rms}}$','$\\gamma_{\mathsf{i,rms}}$','$\\tilde\\gamma_{\mathsf{i,rms}}$','$P_{\mathsf{i,rms}}$']
   
#=================================================================
# Set up figure:
fig = plt.figure(1,figsize=[6.6,6])
ax1 = plt.axes([0.2, 0.2, 0.75, 0.75])

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel(ylabel[field], fontsize=30)

# Loop over directories and plot results:
xmax=0.0
for m,dir in enumerate(dir_list):
   # Select data file:
   dfile='evolution/'+file_prefix[field]+'norms.asc'
   # Open file and read:
   in_file=open(dir+dfile,'r')
   time, y1, y2, y3 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   if option==3:
      # Remove initial time if showing imbalanced data:
      time=np.delete(time,0)
      y1=np.delete(y1,0)
      y2=np.delete(y2,0)
      y3=np.delete(y3,0)
   xmax=max(xmax,time[-1])
   y=np.array([y1,y2,y3])
   
   # Plot:
   if leg==1:
      ax1.plot(time,y[option-1],c=colorlist[m],lw=2,label=label_list[m])
      ax1.legend(loc='upper right',prop={'size':20})
   else:
      ax1.plot(time,y[option-1],c=colorlist[m],lw=2)

#=========================================================================
# Determine nice x & y ranges for plotting:
ax1.set_xlim(0.0,xmax)
ax1.set_yscale('log')
ax1.set_box_aspect(1)

# Save figure:
fig.savefig(outfile, format='eps', dpi=1200)
#fig.savefig(outfile, format='eps', bbox_inches='tight', dpi=1200)

print()
print(' To view the image, type')
print()
print(' gv',outfile,'&')
print()
