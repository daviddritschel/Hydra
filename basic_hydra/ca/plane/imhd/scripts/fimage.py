#!/usr/bin/env python

#==========================================================================
#   Imaging routine superposing field lines (contours of the potential A)
#   on the (potential) vorticity q at a selected time frame.

#   Also creates an image of the current density j.

#   Reads data from qq.r4, jj.r4 and aa.r4.  The latter must be first 
#   generated using magpot.

#   Writes qnnnn.png & jnnnn.png where nnnn = 0000, 0001 ... is the time
#   frame (note 0000 corresponds to t = 0 or fortran record 1).

#   Adapted from ~/clouds/mpic/cimage.py by D G Dritschel near Newcastle
#   on 5/9/2017.
#==========================================================================

import os
import numpy as np

lplotting_cb=True
try:
   from plotting_cb import *
except:
   lplotting_cb=False

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib as mpl
import matplotlib.colors as clrs
import matplotlib.cm as cm
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

#--------------------------------------------------------------
## Global settings

# Set tick label size:
label_size = 30
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# Set x tick width and size:
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1
# Set y tick width and size:
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 1
# Ensure negative contour levels are solid as well:
mpl.rcParams['contour.negative_linestyle'] = 'solid'

#--------------------------------------------------------------------
# Define various useful functions:

def read_file(filename,nx,ny,frame):
   # Read in data from the file in filename at a selected time frame:
   in_file = open(filename,'r')
   raw_array = np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   ngp=nx*ny
   return raw_array[frame*(ngp+1)+1:(frame+1)*(ngp+1)].reshape(ny,nx).T

def contint(fmin,fmax):
   # Determines a nice contour interval (giving 10-20 divisions with
   # interval 1, 2 or 5x10^m for some m) given the minimum & maximum
   # values of the field data (fmin & fmax).

   fmax=0.9999999*fmax
   fmin=0.9999999*fmin
   # The 0.99... factor avoids having a superfluous tick interval
   # in cases where fmax-fmin is 10^m or 2x10^m

   mpow=0
   rmult=fmax-fmin
   while rmult < 10.0:
      mpow+=1
      rmult=rmult*10.0

   while rmult >= 100.0:
      mpow-=1
      rmult=rmult/10.0

   emag=10.0**(float(-mpow))

   kmult=int(rmult/10.0)

   if kmult < 1:
      ci=emag
   elif kmult < 2:
      ci=2.0*emag
   elif kmult < 4:
      ci=4.0*emag
   elif kmult < 8:
      ci=10.0*emag
   else:
      ci=20.0*emag

   return ci

#===========================================================================
#                       Main programme starts here
#===========================================================================
if __name__ == "__main__":

   # Get image dimensions, data save interval and mean magnetic field
   # from job_info file:
   f_info = open('job_info','r')

   for j in range(8):
      line=f_info.readline()
   string=f_info.readline().split()
   nx=int(string[-1])
   ny=nx

   for j in range(3):
      line=f_info.readline()
   string=f_info.readline().split()
   dtsave=float(string[-1])

   for j in range(16):
      line=f_info.readline()
   string=f_info.readline().split()
   B0=float(string[-1])

   f_info.close()

   kfr_in = input('Time frame to image (0 for first; default 20)? ')
   kfr = int(kfr_in or 20)

   time=dtsave*float(kfr)
   print()
   print(' *** This corresponds to t = ',time)
   print()

   # Automatically generate output filename:
   if kfr < 10:
      frame='000'+str(kfr)
   elif kfr < 100:
      frame='00'+str(kfr)
   elif kfr < 1000:
      frame='0'+str(kfr)
   else:
      frame=str(kfr)

   file1='q'+frame+'.png'
   file2='j'+frame+'.png'

   # Read (potential) vorticity:
   Q=read_file('evolution/qq.r4',nx,ny,kfr)

   # Read current density:
   C=read_file('evolution/jj.r4',nx,ny,kfr)

   # Read magnetic potential:
   A=read_file('evolution/aa.r4',nx,ny,kfr)

   # Add mean part (B0*y):
   pi=3.14159265358979
   yg=np.linspace(-pi,pi,ny+1)
   for ix in range(nx):
      A[0:ny-1,ix]+=B0*yg[0:ny-1]

   # Select maximum values to display:
   qmax_def=np.amax(abs(Q))
   qmax_in = input('Maximum |q| to show (default '+str(qmax_def)+')? ')
   qmax = float(qmax_in or qmax_def)
 
   amax_def=np.amax(abs(A))
   amax_in = input('Maximum |A| to show (default '+str(amax_def)+')? ')
   amax = float(amax_in or amax_def)

   cmax_def=np.amax(abs(C))
   cmax_in = input('Maximum |j| to show (default '+str(cmax_def)+')? ')
   cmax = float(cmax_in or cmax_def)
 
   # Set the colourmap:
   colourmap='seismic'
   ccmap=cm.get_cmap(colourmap)
        
   opt_in = input('Add a colourbar (default y)? ')
   cbopt = str(opt_in or 'y')

   opt_in = input('Label x & y axes (default n)? ')
   axlabopt = str(opt_in or 'n')

   aspect=1

   #------------------------------------------------------------------
   # Set up contour levels:
   dq_def=2.0*contint(0.0,qmax)
   dq_in = input('q contour interval (default '+str(dq_def)+')? ')
   dq = float(dq_in or dq_def)
   qmax=dq*float(int( qmax/dq))
   qmin=-qmax
   qncint=int((qmax-qmin)/dq)
   qclevels=np.linspace(qmin,qmax,qncint+1)

   da_def=2.0*contint(0.0,amax)
   da_in = input('A contour interval (default '+str(da_def)+')? ')
   da = float(da_in or da_def)
   amax=da*(float(int( amax/da))+0.5)
   amin=-amax
   ancint=int((amax-amin)/da)
   aclevels=np.linspace(amin,amax,ancint+1)

   dc_def=2.0*contint(0.0,cmax)
   dc_in = input('j contour interval (default '+str(dc_def)+')? ')
   dc = float(dc_in or dc_def)
   cmax=dc*float(int( cmax/dc))
   cmin=-cmax
   cncint=int((cmax-cmin)/dc)
   cclevels=np.linspace(cmin,cmax,cncint+1)

   #------------------------------------------------------------------
   # Plot q and A:
   fig1 = plt.figure(1,figsize=[10,10])
   ax1 = fig1.add_subplot(111)

   ax1.set_xlim(-pi,pi)
   ax1.set_ylim(-pi,pi)
   if axlabopt=='y':
      ax1.set_xlabel('$x$', weight='bold', fontsize=40)
      ax1.set_ylabel('$y$', weight='bold', fontsize=40)
      ax1.set_xticks([-3,-2,-1,0,1,2,3])
      ax1.set_yticks([-3,-2,-1,0,1,2,3])
   else:
      ax1.axes.get_xaxis().set_ticks([])
      ax1.axes.get_yaxis().set_ticks([])

   # Draw image of potential vorticity, q:
   im1=ax1.imshow(Q,cmap=ccmap,vmin=qmin,vmax=qmax,aspect=aspect,extent=(-pi,pi,-pi,pi),origin='lower',interpolation='bilinear')

   # Superpose A contours (field lines):
   ax1.contour(A,aclevels,colors='k',origin='lower',extent=(-pi,pi,-pi,pi),linewidths=2)

   # Optionally add a colourbar:
   if cbopt=='y':
      divider = make_axes_locatable(ax1)
      cax = divider.append_axes("right", size="4%", pad=0.1)
      cbar=fig1.colorbar(im1, cax=cax, ticks=qclevels)
      setp(cbar.ax.yaxis.set_ticklabels(qclevels), weight='bold', fontsize=24)

   # Save figure after cropping:
   fig1.savefig(file1, bbox_inches='tight', pad_inches = 0.025, dpi=100)    

   #------------------------------------------------------------------
   # Plot j:
   fig2 = plt.figure(2,figsize=[10,10])
   ax2 = fig2.add_subplot(111)

   ax2.set_xlim(-pi,pi)
   ax2.set_ylim(-pi,pi)
   if axlabopt=='y':
      ax2.set_xlabel('$x$', weight='bold', fontsize=40)
      ax2.set_ylabel('$y$', weight='bold', fontsize=40)
      ax2.set_xticks([-3,-2,-1,0,1,2,3])
      ax2.set_yticks([-3,-2,-1,0,1,2,3])
   else:
      ax2.axes.get_xaxis().set_ticks([])
      ax2.axes.get_yaxis().set_ticks([])

   # Draw image of current density, j:
   im2=ax2.imshow(C,cmap=ccmap,vmin=cmin,vmax=cmax,aspect=aspect,extent=(-pi,pi,-pi,pi),origin='lower',interpolation='bilinear')

   # Optionally add a colourbar:
   if cbopt=='y':
      divider = make_axes_locatable(ax2)
      cax = divider.append_axes("right", size="4%", pad=0.1)
      cbar=fig2.colorbar(im2, cax=cax, ticks=cclevels)
      setp(cbar.ax.yaxis.set_ticklabels(cclevels), weight='bold', fontsize=24)

   #------------------------------------------------------------------
   # Save figure after cropping:
   fig2.savefig(file2, bbox_inches='tight', pad_inches = 0.025, dpi=100)    

   print()
   print(' *** Created '+file1+' and '+file2)
   print()
