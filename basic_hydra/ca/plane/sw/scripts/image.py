#!/usr/bin/env python

#======================================================================
#   Imaging routine for data produced by sw caps codes

#   Reads data in the current directory
#======================================================================

import os,warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
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

# ==================  Function Definitions ==================

def get_ng():
    # This routine reads the src/parameters.f90 file
    # to obtain the grid resolution of the run ng.
    # This is done by using a simple regexp.
    # It returns an integer with value ng.

    import re

    # Open input file:
    try:
        in_file = open('src/parameters.f90','r')# try opening filename  
    except IOError, message:# error if file not found 
        print >> sys.stderr, ' File could not be opened', message
        sys.exit()

    # Read the file into a list of strings and close it:
    param_file=in_file.readlines()
    in_file.close()

    # Make a list of all strings matching the given regexp: 
    comprehend=[]
    for line in param_file:
        comprehend.append(re.findall('(?<=ng=)\d+',line))
    # Trim the list, select an element, convert to integer and return it:
    comprehend=[comp for comp in comprehend if comp!=[]]
    return int(comprehend[0][0])

def get_cof():
    # This routine reads the src/parameters.f90 file
    # to obtain the Coriolis frequency f.
    # This is done by using a simple regexp.
    # It returns a floating point number with value f.

    import re

    # Open input file:
    try:
        in_file = open('src/parameters.f90','r')# try opening filename  
    except IOError, message:# error if file not found 
        print >> sys.stderr, ' File could not be opened', message
        sys.exit()

    # Read the file into a list of strings and close it:
    param_file=in_file.readlines()
    in_file.close()

    # Make a list of all strings matching the given regexp: 
    comprehend=[]
    for line in param_file:
        comprehend.append(re.findall('(?<=cof=)\d+',line))
    # Trim the list, select an element, convert to integer and return it:
    comprehend=[comp for comp in comprehend if comp!=[]]
    return float(comprehend[0][0])

def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

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
    
def make_plot(Z,cbopt):
    #------------------------------------------------------------------
    # Get contour intervals:
    zmax_def=abs(Z).max()
    zmax=float(raw_input('Maximum magnitude of field to show? (default '+str(zmax_def)+') ') or zmax_def)
    zmin=-zmax

    dZ=contint(zmin,zmax)
    jmax=int(zmax/dZ)
    clevels=np.linspace(-dZ*float(jmax),dZ*float(jmax),2*jmax+1)
    print ' Contour levels:'
    print clevels

    #------------------------------------------------------------------
    # Set up figure:

    # define a square figure:
    fig1 = plt.figure(1,figsize=[8,8])
    ax1 = fig1.add_subplot(111)
    
    #------------------------------------------------------------------
    # Add axes and tick labels:
#    ax1.set_xlabel('$x$', fontsize=25)
    ax1.set_xticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
    ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=24)
#    ax1.set_ylabel('$y$', fontsize=25)
    ax1.set_yticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
    ax1.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=24)

    # Draw image:
    im1=ax1.imshow(Z,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')

    # Superpose contours:
    #ax1.contour(Z,clevels,colors='k',origin='lower',extent=(0.0,1.0,0.0,1.0),linewidths=2)

    # Optionally add a colourbar:
    if cbopt=='y':
       divider = make_axes_locatable(ax1)
       cax = divider.append_axes("right", size="4%", pad=0.1)
       cbar=fig1.colorbar(im1, cax=cax, ticks=clevels)
       #setp(cbar.ax.yaxis.get_ticklabels(), weight='bold', fontsize=24)  
       setp(cbar.ax.yaxis.set_ticklabels(clevels), weight='bold', fontsize=24)

    #------------------------------------------------------------------
    # Save figure after cropping:
    plt.savefig(filename+'.eps', format='eps', dpi=300)    
    print
    print ' *** To view the image, type'
    print
    print ' gv '+filename+'.eps'
    print
    
if __name__ == "__main__":
    #==========================================================================
    # Select field to display:
    field=str(raw_input('Field to image (dd, gg, hh, qq or zz)? (default zz) ') or 'zz')

    if field == 'dd':
        #Get Coriolis frequency:
        cof=get_cof()
        scf=1.0/cof
        print 'Note: the divergence is scaled by the Coriolis frequency.'
    elif field == 'gg':
        #Get Coriolis frequency:
        cof=get_cof()
        scf=1.0/cof**2
        print 'Note: the acceleration divergence is scaled by the Coriolis frequency squared.'
    elif field == 'hh':
        scf=1.0
        print 'Note: hh is the dimensionless height anomaly.'
    elif field == 'qq':
        #Get Coriolis frequency:
        cof=get_cof()
        scf=1.0/cof
        print 'Note: the PV is scaled by the Coriolis frequency.'
    elif field == 'zz':
        #Get Coriolis frequency:
        cof=get_cof()
        scf=1.0/cof
        print 'Note: the relative vorticity is scaled by the Coriolis frequency.'
    else:
        print ' *** Not an available field *** EXITING!'
        quit()

    kfr=int(raw_input('Time frame to image? (default 100) ') or 100)
    if kfr < 10:
       frame='000'+str(kfr)
    elif kfr < 100:
       frame='00'+str(kfr)
    elif kfr < 1000:
       frame='0'+str(kfr)
    else:
       frame=str(kfr)
    
    #Output filename:
    filename=field+frame

    #Get resolution:
    ng=get_ng()

    #Read selected frame (kfr) from the data:
    in_file=open(field+'.r4','r')
    raw_array=np.fromfile(in_file,dtype=np.float32)
    in_file.close()
    N=ng*ng

    Z=np.empty((ng+1,ng+1))
    Z[0:ng,0:ng]=scf*raw_array[kfr*(N+1)+1:(kfr+1)*(N+1)].reshape(ng,ng).T

    # Add periodic edges:
    Z[ng,0:ng]=Z[0,0:ng]
    Z[0:ng+1,ng]=Z[0:ng+1,0]

    cbopt=str(raw_input('Add a colourbar? (default y) ') or 'y')

    make_plot(Z,cbopt)
