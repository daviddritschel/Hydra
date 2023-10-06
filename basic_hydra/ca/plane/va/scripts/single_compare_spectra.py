#!/usr/bin/env python

# This script plots spectra for either height, vorticity, divergence
# or acceleration divergence from data in directories specified below.

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
# Corresponding labels on the curves plotted (modify as necessary):
label_list=['$H=0$','$H=0.1$','$H=0.2$','$H=0.4$']
# Corresponding colours (allow up to 4 directories):
colorlist=['k','b','r','m']

#=================================================================
# Select data to plot:
print()
print(' This script compares spectra of different types in different directories.')
print()
print(' Show spectra for')
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

opt_in = input(' show (1) full, (2) balanced or (3) imbalanced spectra (default 1)? ')
option = int(opt_in or 1)

print()
t_in = input(' Time to show (default 500)? ')
t = float(t_in or 500.0)

print()
opt_in = input(' Add a legend (1 = yes, 0 = no; default 0)? ')
leg = int(opt_in or 0)
   
file_prefix=['h','z','d','g','t','p']

if option==1:
   bal=''
   file_suffix='spec.asc'
   ylabel=['$\log_{10}S_{h}$','$\log_{10}S_{\\zeta}$','$\log_{10}S_{\\delta}$','$\log_{10}S_{\\gamma}$','$\log_{10}S_{\\tilde\\gamma}$','$\log_{10}S_{P}$']
elif option==2:
   bal='b'
   file_suffix='bspec.asc'
   ylabel=['$\log_{10}S_{h_{\mathsf{b}}}$','$\log_{10}S_{\\zeta_{\mathsf{b}}}$','$\log_{10}S_{\\delta_{\mathsf{b}}}$','$\log_{10}S_{\\gamma_{\mathsf{b}}}$','$\log_{10}S_{\\tilde\\gamma_{\mathsf{b}}}$','$\log_{10}S_{P_{\mathsf{b}}}$']
else:
   bal='i'
   file_suffix='ispec.asc'
   ylabel=['$\log_{10}S_{h_{\mathsf{i}}}$','$\log_{10}S_{\\zeta_{\mathsf{i}}}$','$\log_{10}S_{\\delta_{\mathsf{i}}}$','$\log_{10}S_{\\gamma_{\mathsf{i}}}$','$\log_{10}S_{\\tilde\\gamma_{\mathsf{i}}}$','$\log_{10}S_{P_{\mathsf{i}}}$']

outfile=bal+file_prefix[field]+'_spectra_t'+str(int(t+0.01))+'.eps'

#=================================================================
# Set up figure:
fig = plt.figure(1,figsize=[6.6,6])
ax1 = plt.axes([0.2, 0.2, 0.75, 0.75])

ax1.set_xlabel('$\log_{10}k$', fontsize=30)
ax1.set_ylabel( ylabel[field], fontsize=30)

# Loop over directories and plot results:
smin=32.0
smax=-32.0
xmax=0.0
for m,dir in enumerate(dir_list):
   # Select data file:
   dfile='spectra/'+file_prefix[field]+file_suffix
   # Open file and read:
   in_file=open(dir+dfile,'r')
   first_line=in_file.readline()
   nx=int(first_line.split()[-1])
   in_file.seek(0)

   kc=int(2.0*float(nx)/3.0)
   xmax=max(xmax,0.05*int(20.0*np.log10(float(kc))+1.0))

   # Read in the full data to a 1d array and close input file:
   raw_data = np.fromfile(file=in_file,dtype=float,sep='\n')
   in_file.close()

   # Determine the number of frames:
   nframes = int(len(raw_data)/(2*nx+2))  

   # Shape the data array into a useful shape for plotting:
   frames=range(0,nframes)
   time=[raw_data[i*(2*nx+2)] for i in frames]
   dt=time[1]-time[0]
   tim_eles = [i*(2*nx+2)+j for i in frames for j in range(2)]
   shaped_data = np.delete(raw_data,tim_eles)[0:(2*nx+2)*nframes].reshape((nframes,nx,2))
   k=np.zeros((nframes,nx))
   h=np.zeros((nframes,nx))
   for i in frames:
      k[i,:]=shaped_data[i].transpose()[0][0:nx]
      h[i,:]=shaped_data[i].transpose()[1][0:nx]

   # Select frame:
   ic=int(t/dt+0.01)

   # Find min & max values:
   smin=min(smin,np.amin(h[ic,:kc]))
   smax=max(smax,np.amax(h[ic,:kc]))

   # Plot:
   if leg==1:
      ax1.plot(k[ic,:kc],h[ic,:kc],c=colorlist[m],lw=2,label=label_list[m])
      ax1.legend(loc='lower center',prop={'size':20})
   else:
      ax1.plot(k[ic,:kc],h[ic,:kc],c=colorlist[m],lw=2)
   
#=========================================================================
# Determine nice x & y ranges for plotting:
ax1.set_xlim(0.0,xmax)

ymin=0.5*int(2.0*(smin+1000.0))-1000.0
ymax=0.5*int(2.0*(smax+1000.0)+1.0)-1000.0
ax1.set_ylim(ymin,ymax)
ax1.set_aspect(xmax/(ymax-ymin))

# Save figure:
fig.savefig(outfile, format='eps', dpi=1200)
#fig.savefig(outfile, format='eps', bbox_inches='tight', dpi=1200)

print()
print(' To view the image, type')
print()
print(' gv',outfile)
print()
