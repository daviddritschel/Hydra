#!/usr/bin/env python

# This plots the contours in pvcont.dat, e.g. those previously generated
# by "tracer" (see post/tracer.f90).

# Written 4 November 2017 by D G Dritschel @ St Andrews

#=====perform the various imports========
import warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rcParams

## global plot settings
rcParams.update({'figure.autolayout': True})

# set tick label size:
label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 8
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1

# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 8
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 1

warnings.simplefilter("ignore",DeprecationWarning)

#==========================================================
# Set default plot parameters:

# Maximum colour saturation level (1 = black):
satmax=0.9

# Figure width in inches:
figwidth=8.0

# Line width of plotted contours (integer):
linewidth=2

# resolution of saved figure (dpi, integer):
figres=200

#==========================================================
# Open input data file:
in_file=open('pvcont.dat','r')

nfr_in=input('Frame to plot (1 for the first - default 1)? ')
nfr=int(nfr_in or 1)
xmax_in=input('Maximum |x|,|y| to show (default 1.0)? ')
xmax=float(xmax_in or 1.0)
xmin=-xmax
ymax=xmax
ymin=-ymax
ntint_in=input('Number of tick intervals on each axis (default 4)? ')
ntint=int(ntint_in or 4)

if nfr > 1:
   # Skip frames:
   for k in range(1,nfr):
      record=in_file.readline()
      nc=int(record.split()[0])
      npt=int(record.split()[1])
      for i in range(1,nc+npt+1):
         record=in_file.readline()

# Read the chosen frame and plot:
record=in_file.readline()
nc=int(record.split()[0])
npt=int(record.split()[1])
t=float(record.split()[2])

print()
print(' nc = ',nc,'   npt = ',npt,'   t = ',t)

#========================================================
# Create figure:
fig = plt.figure(1,figsize=[figwidth,figwidth*ymax/xmax])
ax = plt.axes([0.2, 0.2, 0.7, 0.7])

# Fix domain of view:
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])

# Plot ticks on both axes in each direction:
ax.tick_params(top=True, right=True)

# Force chosen tick interval on each axis:
dxtick=(xmax-xmin)/float(ntint)
ax.xaxis.set_ticks(np.arange(xmin, xmax+dxtick, dxtick))

dytick=(ymax-ymin)/float(ntint)
ax.yaxis.set_ticks(np.arange(ymin, ymax+dytick, dytick))

# Loop over contours and plot:
np=[]
ind=[]
# Get minimum and maximum "levels" for choosing contour colours:
indmin=1000000
indmax=-1000000
for j in range(0,nc):
   record=in_file.readline()
   np.append(int(record.split()[0]))
   lev=int(record.split()[2])
   ind.append(lev)
   indmin=min(indmin,lev)
   indmax=max(indmax,lev)

for j in range(0,nc):
   if indmin == indmax:
      # Only one level is present: use a black line:
      color='k'
   else:
      # Choose contour colour depending on sign and magnitude of ind:
      if ind[j] > 0:
         color=plt.cm.Reds(satmax*float(ind[j])/float(indmax))
      elif ind[j] < 0:
         color=plt.cm.Blues(satmax*float(ind[j])/float(indmin))

   x=[]
   y=[]

   record=in_file.readline()
   x0=float(record.split()[0])
   y0=float(record.split()[1])
   x.append(x0)
   y.append(y0)

   for i in range(1,np[j]):
      record=in_file.readline()
      x.append(float(record.split()[0]))
      y.append(float(record.split()[1]))

   x.append(x0)
   y.append(y0)

   ax.plot(x,y,color=color,lw=linewidth)

# Save figure after cropping:
plt.savefig('contours.png', bbox_inches='tight', pad_inches = 0.025, dpi=figres)

plt.show()

print()
print(' *** Image saved in contours.png')
print()
