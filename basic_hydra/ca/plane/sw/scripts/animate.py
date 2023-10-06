#!/usr/bin/env python

#=====================================================================
#   Animates data in any .r8 file over a selected time range

#   Written by D G Dritschel, 27 November 2018 @ St Andrews
#=====================================================================

#========== Perform the generic imports =========
import sys,os,warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.animation as anim
import matplotlib.colors as clrs
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

# set tick label size:
label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
# set axes width:
mpl.rcParams['axes.linewidth'] = 2
#=========================================

#-----------------------------------------------------------
# Initialisation
#-----------------------------------------------------------
# Figure width in inches:
width=6.0

# Get image dimensions from qq_init.r8:
myx=(os.path.getsize('qq_init.r8')-8)/8
nx=int(np.sqrt(float(myx)+1.e-6))
ny=nx

length=nx*ny+1
frame_bytes=8*length

prefix=str(raw_input('File prefix (before .r8, default qq)? ') or 'qq')
dataset=prefix+'.r8'

# Get final frame from chosen dataset:
record2=os.path.getsize(dataset)/frame_bytes-1

record1=int(raw_input('First record to read (default 0)? ') or 0)
record2=int(raw_input(' Last record to read (default '+str(record2)+')? ') or record2)

# Set up figure:
fig = plt.figure(1,figsize=[width,width])
ax = fig.add_subplot(111)

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

# Formatting for displayed time:
#time_template = 't = %.1f'
time_template = 't = %.0f'

# Set up list of images for animation:
ims=[]

#-----------------------------------------------------------
# Read selected range of frames and process:
raw_array=np.empty(length, dtype=np.float64)
with open(dataset,'r') as file:

   for record in range(record1,record2+1):
      file.seek(frame_bytes * record)
      bytes=file.read(frame_bytes)
      raw_array=np.frombuffer(bytes, dtype=np.float64).copy()

      t=float(raw_array[0])
      print ' t = ',t

      Z=np.empty((nx+1,ny+1))
      Z[0:nx,0:ny]=raw_array[1:nx*ny+1].reshape(nx,ny)
      Z[nx,0:ny]=Z[0,0:ny]
      Z[0:nx+1,ny]=Z[0:nx+1,0]

      zmin=np.amin(Z)
      zmax=np.amax(Z)

      # Plot frame:
      im1 = ax.imshow(Z.T,cmap=cm.terrain,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
      # Add time as a text string:
      im2 = ax.text(0.025, 0.955, time_template%(t), fontsize=16, weight='bold', transform=ax.transAxes)
      # Add to list of images:
      ims.append([im1,im2])

#-----------------------------------------------------------
# Run animation:
ani = anim.ArtistAnimation(fig, ims, interval=200, repeat_delay=2000, blit=False)

movie_file=prefix+str(record1)+'-'+str(record2)+'.mp4'
print
print ' Creating movie ('+movie_file+') ... (this can take some time)'
ani.save(movie_file, fps=10, extra_args=['-vcodec', 'libx264'])

print
print ' Display the movie by typing'
print
print ' mplayer '+movie_file+' -loop 0'
print
