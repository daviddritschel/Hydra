#!/usr/bin/env python

# This script plots the buoyancy and vorticity spectra, comparing 
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
print()
t_in = input(' Time to show (default 900)? ')
t = float(t_in or 900)

# Open ene.asc file in one directory to get time between frames:
in_file=open(dir_list[0]+'/evolution/ecomp.asc','r')
time, ekin, epot, etot = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
nt=len(time)

dt=time[-1]/float(nt-1)
# Frame corresponding to time chosen:
frame=int(t/dt+0.5)
print()
print(' Plotting spectra at time',time[frame])

# Set up figure:
fig, (ax1, ax2) = plt.subplots(figsize=[12.5,6], nrows=1, ncols=2)

ax1.set_xlabel('$\\log_{10}k$', fontsize=30)
ax1.set_ylabel('$\\log_{10}S_b$', fontsize=30)
#ax1.yaxis.set(major_locator=MultipleLocator(0.1), \
#                 major_formatter=FormatStrFormatter('%1.1f'))
ax1.set_box_aspect(1)

ax2.set_xlabel('$\\log_{10}k$', fontsize=30)
ax2.set_ylabel('$\\log_{10}S_{\\zeta}$', fontsize=30)
#ax2.yaxis.set(major_locator=MultipleLocator(0.1), \
#                 major_formatter=FormatStrFormatter('%1.1f'))
ax2.set_box_aspect(1)

#=================================================================
# Loop over directories and plot results:
s1max=-1000.0
s2max=-1000.0
for m,dir in enumerate(dir_list):
   # Read and plot b spectrum:
   in_file1=open(dir+'/spectra/bspec.asc','r')
   first_line=in_file1.readline()
   kmax=int(first_line.split()[-1])
   in_file1.seek(0)
   nx=kmax

   raw_data1 = np.fromfile(file=in_file1,dtype=float,sep='\n')
   in_file1.close()
   nframes = int(len(raw_data1)/(2*kmax+2))  

   frames=range(0,nframes)
   time=[raw_data1[i*(2*kmax+2)] for i in frames]
   tim_eles = [i*(2*kmax+2)+j for i in frames for j in range(2)]
   shaped_data = np.delete(raw_data1,tim_eles)[0:(2*kmax+2)*nframes].reshape((nframes,kmax,2))
   k=shaped_data[frame].transpose()[0][0:nx]
   s1=shaped_data[frame].transpose()[1][0:nx]
   s1max=max(s1max,np.amax(s1))

   ax1.plot(k,s1,c=colorlist[m],lw=2,label=label_list[m])
   ax1.legend(loc='lower left',prop={'size':20})
   
   # Read and plot zeta spectrum:
   in_file2=open(dir+'/spectra/zspec.asc','r')
   first_line=in_file2.readline()
   kmax=int(first_line.split()[-1])
   in_file2.seek(0)
   nx=kmax

   raw_data2 = np.fromfile(file=in_file2,dtype=float,sep='\n')
   in_file2.close()
   nframes = int(len(raw_data2)/(2*kmax+2))  

   frames=range(0,nframes)
   time=[raw_data2[i*(2*kmax+2)] for i in frames]
   tim_eles = [i*(2*kmax+2)+j for i in frames for j in range(2)]
   shaped_data = np.delete(raw_data2,tim_eles)[0:(2*kmax+2)*nframes].reshape((nframes,kmax,2))
   k=shaped_data[frame].transpose()[0][0:nx]
   s2=shaped_data[frame].transpose()[1][0:nx]
   s2max=max(s2max,np.amax(s2))

   ax2.plot(k,s2,c=colorlist[m],lw=2)
   
s1max=0.2*int(5.0*(s1max+1000.0)+1.0)-1000.0
s2max=0.2*int(5.0*(s2max+1000.0)+1.0)-1000.0

ax1.set_ylim(s1max-4.0,s1max)
ax2.set_ylim(s2max-2.0,s2max)

#=========================================================================
# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0)

fname='spectra'+str(t)+'.eps'
fig.savefig(fname, format='eps', dpi=1200)

print()
print(' To view the image, type')
print()
print(' gv '+fname+' &')
print()
