#!/usr/bin/env python

# This script plots qbar & q_e vs phi as well as ubar vs phi side by side

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
fig, (ax1, ax2) = plt.subplots(figsize=[12.1,6.0], nrows=1, ncols=2)

ax1.set_xlabel('$\\bar{q}$, $q_e$', fontsize=30)
ax1.set_ylabel('$\phi$', fontsize=30)
ax1.set_ylim(-np.pi/2.0,np.pi/2.0)
ax1.axvline(0.0,color='b',ls='--',lw=1)
ax1.yaxis.set_ticks([-np.pi/2.0,-np.pi/4.0,0.0,np.pi/4.0,np.pi/2.0])
ax1.set_yticklabels([r'$-\pi/2$',r'$-\pi/4$',r'$0$',r'$\pi/4$',r'$\pi/2$'],fontsize=20)

ax2.set_xlabel('$\\bar{u}$', fontsize=30)
#ax2.set_ylabel('$\phi$', fontsize=30)
ax2.set_ylim(-np.pi/2.0,np.pi/2.0)
ax2.axvline(0.0,color='b',ls='--',lw=1)
ax2.yaxis.set_ticks([-np.pi/2.0,-np.pi/4.0,0.0,np.pi/4.0,np.pi/2.0])
ax2.set_yticklabels([r'$-\pi/2$',r'$-\pi/4$',r'$0$',r'$\pi/4$',r'$\pi/2$'],fontsize=20)

#=================================================================
# Open file containing qbar:
in_file=open('evolution/zq.asc','r')
# Read the first header line to get ng:
first_line=in_file.readline()
ng=int(first_line.split()[-1])
in_file.seek(0)

# Read in the full data to a 1d array and close input file:
raw_data = np.fromfile(file=in_file,dtype=float,sep='\n')
in_file.close()
nframes = int(len(raw_data)/(2*ng+2))  

# Shape the data array for plotting:
frames=range(0,nframes)
time=[raw_data[i*(2*ng+2)] for i in frames]
tim_eles = [i*(2*ng+2)+j for i in frames for j in range(2)]
shaped_data = np.delete(raw_data,tim_eles)[0:(2*ng+2)*nframes].reshape((nframes,ng,2))
x1=np.zeros((nframes,ng))
y1=np.zeros((nframes,ng))
for i in frames:
   x1[i,:]=shaped_data[i].transpose()[0][0:ng]
   y1[i,:]=shaped_data[i].transpose()[1][0:ng]

#--------------------------------------------------------
# Open file containing ubar:
in_file=open('evolution/zu.asc','r')

# Read in the full data to a 1d array and close input file:
raw_data = np.fromfile(file=in_file,dtype=float,sep='\n')
in_file.close()

# Shape the data array for plotting:
shaped_data = np.delete(raw_data,tim_eles)[0:(2*ng+2)*nframes].reshape((nframes,ng,2))
x3=np.zeros((nframes,ng))
y3=np.zeros((nframes,ng))
for i in frames:
   x3[i,:]=shaped_data[i].transpose()[0][0:ng]
   y3[i,:]=shaped_data[i].transpose()[1][0:ng]

# Remove solid body rotation in u_bar (in the array x3):
f1112=1.0/12.0
clat=np.cos(y3[0,:])
rsum=f1112*(clat[0]+clat[ng-1])+np.sum(clat[1:ng-1])
rsum=1.0/rsum
for i in frames:
   ombar=(f1112*(x3[i,0]+x3[i,ng-1])+np.sum(x3[i,1:ng-1]))*rsumi
   x3[i,:]=x3[i,:]-ombar*clat
   
#--------------------------------------------------------
# Open file containing q_e:
in_file=open('evolution/qe.asc','r')

# Read in the full data to a 1d array and close input file:
raw_data = np.fromfile(file=in_file,dtype=float,sep='\n')
in_file.close()
ngp1=ng+1
nframes = int(len(raw_data)/(2*ngp1+2))  

# Shape the data array for plotting:
frames=range(0,nframes)
time=[raw_data[i*(2*ngp1+2)] for i in frames]
tim_eles = [i*(2*ngp1+2)+j for i in frames for j in range(2)]
shaped_data = np.delete(raw_data,tim_eles)[0:(2*ngp1+2)*nframes].reshape((nframes,ngp1,2))
x2=np.zeros((nframes,ngp1))
y2=np.zeros((nframes,ngp1))
for i in frames:
   x2[i,:]=shaped_data[i].transpose()[0][0:ngp1]
   y2[i,:]=shaped_data[i].transpose()[1][0:ngp1]

#--------------------------------------------------------
tmax=time[-1]
print()
t_in = input(' Time to show (default '+str(tmax)+')? ')
t = int(t_in or tmax)
ic=int(t/(time[1]-time[0])+0.00001)
t=time[ic]

ax1.plot(4.0*np.pi*np.sin(y1[ic]),y1[ic],c='m',lw=1,label='$2\\Omega\\sin\\phi$')
ax1.plot(x1[ic],y1[ic],c='k',lw=1,label='$\\bar{q}$')
#ax1.plot(x2[ic],y2[ic],c=(0.4,0.4,0.4),lw=1,label='$q_e$')
ax1.plot(x2[ic],y2[ic],c='r',lw=1,label='$q_e$')
ax2.plot(x3[ic],y3[ic],c='k',lw=1)

ax1.legend(loc='upper left',prop={'size':20})

#plt.setp(ax2.get_yticklabels(), visible=False)

# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0.1)

ofile='stair_t'+str(t)+'.eps'
fig.savefig(ofile, format='eps', dpi=300)

print(' To display the results, type:')
print()
print(' gv',ofile)
print()
