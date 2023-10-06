#!/usr/bin/env python

# This script plots spectra for height, vorticity, divergence and 
# acceleration divergence from data in spectra.asc and alt-spectra.asc 
# in separate SW directories specified below.

#========== Perform the generic imports =========
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
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
t=float(raw_input('Time to show (default 25)? ') or 25.0)
print

#=================================================================
# Set up figure:

height=15.1

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=[16,height], nrows=2, ncols=2)

#ax1.set_xlabel('$\log_{10}k$', fontsize=30)
ax1.set_xlim(0.0,2.6)

#ax2.set_xlabel('$\log_{10}k$', fontsize=30)
ax2.set_xlim(0.0,2.6)

ax3.set_xlabel('$\log_{10}k$', fontsize=30)
ax3.set_xlim(0.0,2.6)

ax4.set_xlabel('$\log_{10}k$', fontsize=30)
ax4.set_xlim(0.0,2.6)

#=================================================================
outfile='spectra_swca_t'+str(int(t+0.01))+'res.eps'
datafile='spectra.asc'
altdatafile='alt-spectra.asc'
ax1.set_ylabel('$\log_{10}S_{h}$', fontsize=30)
ax2.set_ylabel('$\log_{10}S_{\\zeta}$', fontsize=30)
ax3.set_ylabel('$\log_{10}S_{\\delta}$', fontsize=30)
ax4.set_ylabel('$\log_{10}S_{\\gamma}$', fontsize=30)
print
# Range in log_10 spectra to show:
yrange=12.0
ymax=float(raw_input('Maximum value in log_10   h   spectrum to show (default -3)? ') or -3.0)
ax1.set_ylim(ymax-yrange,ymax)
yrange=9.0
ymax=float(raw_input('Maximum value in log_10  zeta spectrum to show (default  0)? ') or 0.0)
ax2.set_ylim(ymax-yrange,ymax)
ymax=float(raw_input('Maximum value in log_10 delta spectrum to show (default -4)? ') or -4.0)
ax3.set_ylim(ymax-yrange,ymax)
ymax=float(raw_input('Maximum value in log_10 gamma spectrum to show (default  1)? ') or 1.0)
ax4.set_ylim(ymax-yrange,ymax)

#=================================================================
# List of directories to compare (need the final /):
dir_list=['bal_ng256/','bal_ng512/']
# Corresponding labels on the curves plotted:
label_list=['$n=256$','$n=512$']
# Corresponding line styles:
dashlist=[(1,0.0001),(1,0.0001),(8,3),(2,5)]
# Corresponding colours (or shades of grey):
colorlist=[(0.0,0.0,0.0),(0.4,0.4,0.4),(0.6,0.6,0.6)]

#=================================================================
# Loop over directories and plot results:
for m,dir in enumerate(dir_list):
   # Open input file to read h spectrum:
   in_file=open(dir+altdatafile,'r')
   # Read the first header line to get kmax:
   first_line=in_file.readline()
   kmax=int(first_line.split()[-1])
   in_file.seek(0)

   nx=kmax

   # Read in the full data to a 1d array and close input file:
   raw_data = np.fromfile(file=in_file,dtype=float,sep='\n')
   in_file.close()

   # Determine the number of frames:
   nframes = int(len(raw_data)/(4*kmax+2))  

   # Shape the data array into a useful shape for plotting:
   frames=range(0,nframes)
   time=[raw_data[i*(4*kmax+2)] for i in frames]
   dt=time[1]-time[0]
   tim_eles = [i*(4*kmax+2)+j for i in frames for j in range(2)]
   shaped_data = np.delete(raw_data,tim_eles)[0:(4*kmax+2)*nframes].reshape((nframes,kmax,4))
   h=np.zeros((nframes,nx))
   for i in frames:
      h[i,:]=shaped_data[i].transpose()[1][0:nx]

   # Open input file to read remaining spectra:
   in_file=open(dir+datafile,'r')
   # Read the first header line to get kmax:
   first_line=in_file.readline()
   kmax=int(first_line.split()[-1])
   in_file.seek(0)

   nx=kmax

   # Read in the full data to a 1d array and close input file:
   raw_data = np.fromfile(file=in_file,dtype=float,sep='\n')
   in_file.close()

   # Determine the number of frames:
   nframes = int(len(raw_data)/(4*kmax+2))  

   # Shape the data array into a useful shape for plotting:
   frames=range(0,nframes)
   time=[raw_data[i*(4*kmax+2)] for i in frames]
   dt=time[1]-time[0]
   tim_eles = [i*(4*kmax+2)+j for i in frames for j in range(2)]
   shaped_data = np.delete(raw_data,tim_eles)[0:(4*kmax+2)*nframes].reshape((nframes,kmax,4))
   k=np.zeros((nframes,nx))
   z=np.zeros((nframes,nx))
   d=np.zeros((nframes,nx))
   g=np.zeros((nframes,nx))
   for i in frames:
      k[i,:]=shaped_data[i].transpose()[0][0:nx]
      z[i,:]=shaped_data[i].transpose()[1][0:nx]
      d[i,:]=shaped_data[i].transpose()[2][0:nx]
      g[i,:]=shaped_data[i].transpose()[3][0:nx]

   # Select frame:
   ic=int(t/dt+0.01)

   ax1.plot(k[ic],h[ic],dashes=dashlist[m],c=colorlist[m],lw=3,label=label_list[m])
   ax1.legend(loc='lower left',prop={'size':20}, shadow=True)

   ax2.plot(k[ic],z[ic],dashes=dashlist[m],c=colorlist[m],lw=3)

   ax3.plot(k[ic],d[ic],dashes=dashlist[m],c=colorlist[m],lw=3)

   ax4.plot(k[ic],g[ic],dashes=dashlist[m],c=colorlist[m],lw=3)

#=========================================================================
# Add information about the resolution in the first panel:
ax1.text( 0.05, 0.23, 'SW-CA, $t={x:.0f}$'.format(x=t), transform=ax1.transAxes, fontsize=25)

# Share x axes in top and bottom panels:
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

# Add spacing between panels:
plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)

# Save figure:
fig.savefig(outfile, format='eps', dpi=1200)

print
print ' To view the image, type'
print
print ' gv ',outfile
print
