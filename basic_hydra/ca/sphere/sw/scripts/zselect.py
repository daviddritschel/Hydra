#!/usr/bin/env python

# This script plots the zonal average u & h at a selected time

# *** zonal must be run first to generate the data ***

# ==> Run from the current job directory <==

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

ax1.set_xlabel('$\\bar{u}/\\bar{u}_{\\mathsf{rms}}$,~ $\\bar{h}/\\bar{h}_{\\mathsf{rms}}$', fontsize=30)
ax1.set_ylabel('$\phi$', fontsize=30)
ax1.set_ylim(-np.pi/2.0,np.pi/2.0)
ax1.axvline(0.0,color='b',ls='--',lw=1)
ax1.xaxis.set_ticks([-5,-4,-3,-2,-1,0,1,2,3,4,5])
ax1.yaxis.set_ticks([-np.pi/2.0,-np.pi/4.0,0.0,np.pi/4.0,np.pi/2.0])
ax1.set_yticklabels([r'$-\pi/2$',r'$-\pi/4$',r'$0$',r'$\pi/4$',r'$\pi/2$'],fontsize=20)

ax2.set_xlabel('$\\bar{u}/\\bar{u}_{\\mathsf{rms}}$,~ $\\bar{h}/\\bar{h}_{\\mathsf{rms}}$', fontsize=30)
#ax2.set_ylabel('$\phi$', fontsize=30)
ax2.set_ylim(-np.pi/2.0,np.pi/2.0)
ax2.axvline(0.0,color='b',ls='--',lw=1)
ax2.xaxis.set_ticks([-5,-4,-3,-2,-1,0,1,2,3,4,5])
ax2.yaxis.set_ticks([-np.pi/2.0,-np.pi/4.0,0.0,np.pi/4.0,np.pi/2.0])
ax2.set_yticklabels([r'$-\pi/2$',r'$-\pi/4$',r'$0$',r'$\pi/4$',r'$\pi/2$'],fontsize=20)

#=================================================================
# Open file containing hbar:
in_file=open('evolution/zh.asc','r')
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
x2=np.zeros((nframes,ng))
y2=np.zeros((nframes,ng))
for i in frames:
   x2[i,:]=shaped_data[i].transpose()[0][0:ng]
   y2[i,:]=shaped_data[i].transpose()[1][0:ng]

#--------------------------------------------------------
tmax=time[-1]
tsel=tmax/2
print()
t_in = input(' First time to show (default '+str(tsel)+')? ')
t1 = int(t_in or tsel)
ic=int(t1/(time[1]-time[0])+0.5)
t1=time[ic]

# Compute rms values and scale data:
x1rms=np.sqrt(sum(x1[ic]**2)/ng)
x1[ic]=x1[ic]/x1rms

x2rms=np.sqrt(sum(x2[ic]**2)/ng)
x2[ic]=x2[ic]/x2rms

ax1.set_title('$t = {x:.0f}$'.format(x=t1),size=30)
ax1.plot(x1[ic],y1[ic],c='r',lw=2)
ax1.plot(x2[ic],y2[ic],c='k',lw=2)

#--------------------------------------------------------
tsel=tmax
print()
t_in = input(' Second time to show (default '+str(tsel)+')? ')
t2 = int(t_in or tsel)
ic=int(t2/(time[1]-time[0])+0.5)
t2=time[ic]

# Compute rms values and scale data:
x1rms=np.sqrt(sum(x1[ic]**2)/ng)
x1[ic]=x1[ic]/x1rms

x2rms=np.sqrt(sum(x2[ic]**2)/ng)
x2[ic]=x2[ic]/x2rms

ax2.set_title('$t = {x:.0f}$'.format(x=t2),size=30)
ax2.plot(x1[ic],y1[ic],c='r',lw=2)
ax2.plot(x2[ic],y2[ic],c='k',lw=2)
#ax2.plot(x1[ic],y1[ic],c='r',lw=2,label='$\\bar{h}/\\bar{h}_{\\mathsf{rms}}$')
#ax2.plot(x2[ic],y2[ic],c='k',lw=2,label='$\\bar{u}/\\bar{u}_{\\mathsf{rms}}$')
#ax2.legend(loc='center right',prop={'size':25})

#plt.setp(ax2.get_yticklabels(), visible=False)

# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0.1)

ofile='zhu_t'+str(t1)+'_and_t'+str(t2)+'.eps'
fig.savefig(ofile, format='eps', dpi=300)

print(' To display the results, type:')
print()
print(' gv',ofile)
print()
