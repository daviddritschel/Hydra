#!/usr/bin/env python

#================================================================
#  Contours a given field at a given time
#================================================================

#=====perform various generic imports=====
import warnings,os
import numpy as np

from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib as mpl
import matplotlib.colors as clrs
import matplotlib.cm as cm
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})

## global settings

# set tick label size:
label_size = 25
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

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

# --------------------- Function definitions -----------------------
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

#-------------------------------------------------------------------------
#Get image dimensions (nt x ng):
file_bytes=os.path.getsize('hh_init.r8')-8
ng=int(np.sqrt(file_bytes/16))
nt=2*ng

dl=np.pi/float(ng)
lat=np.linspace((dl-np.pi)/2.0,(np.pi-dl)/2.0,ng)
lat=np.insert(lat,0,-np.pi/2.0)
lat=np.append(lat,np.pi/2.0)
lon=np.linspace(-np.pi,np.pi,nt+1)

X, Y = np.meshgrid(lon, lat)

print()
prefix_in=input(' Data file prefix (hh, dd, gg, qq, ...; default hh)? ')
prefix = str(prefix_in or 'hh')
filename='evolution/'+prefix+'.r4'

#Read in data and find min/max values:
in_file = open(filename,'r')
raw_array = np.fromfile(in_file,dtype=np.float32)
in_file.close()

frame_in=input(' Frame to read (0 for first; default 30)? ')
frame = int(frame_in or 30)

N=ng*nt
Z=np.zeros([ng+2,nt+1])
#Grab frame:
Z[1:ng+1,0:nt]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(nt,ng).T
#Add polar values by extrapolation:
Z[0,0:nt]=1.5*Z[1,0:nt]-0.5*Z[2,0:nt]
Z[ng+1,0:nt]=1.5*Z[ng,0:nt]-0.5*Z[ng-1,0:nt]
#Add periodic edge in longitude:
Z[0:ng+2,nt]=Z[0:ng+2,0]

opt_in=input(' Add axis labels and ticks (y/n; default y)? ')
opt = str(opt_in or 'y')

#-------------------------------------------
#Next plot results:
Zmin=np.amin(Z)
Zmax=np.amax(Z)
print(' Min field value = ',Zmin,'  Max field value = ',Zmax)

if prefix=='hh':
   dZ=120/9522.958
   print()
   print(' *** Contour interval used in the paper = ',dZ)
   print()
elif prefix=='qq':
   dZ=2*np.pi*9.3385e-6/7.292e-5
   print()
   print(' *** Contour interval used in the paper = ',dZ)
   print()
elif prefix=='dd':
   dZ=4*np.pi*0.001
   print()
   print(' *** Contour interval used in the paper = ',dZ)
   print()
elif prefix=='gg':
   dZ=(4*np.pi)**2*0.02
   print()
   print(' *** Contour interval used in the paper = ',dZ)
   print()

Zmax=max(Zmax,abs(Zmin))
Zmin=-Zmax
dZ_def=contint(Zmin,Zmax)/2.0
dZ_in=input(' Contour interval? (default '+str(dZ_def)+') ')
dZ = float(dZ_in or dZ_def)
jmax=int(Zmax/dZ+0.5)
Zmax=dZ*(float(jmax)-0.5)
clevels=np.linspace(-Zmax,Zmax,2*jmax)

fig1 = plt.figure(1,figsize=[10,5])
ax1 = fig1.add_subplot(111)

ax1.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])

ax1.yaxis.set_ticks([-np.pi/2.0,0.0,np.pi/2.0])

if opt=='y':
   ax1.set_xlabel('$\\lambda$', fontsize=30)
   ax1.set_ylabel('$\\phi$', fontsize=30)
   ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax1.set_yticklabels([r'$-\pi/2$',r'$0$',r'$\pi/2$'],fontsize=20)
else:
   plt.setp(ax1.get_xticklabels(), visible=False)
   plt.setp(ax1.get_yticklabels(), visible=False)

# Contour:
ax1.contour(X,Y,Z,clevels,colors='k',origin='lower',extent=(-np.pi,np.pi,-np.pi/2,np.pi/2),linewidths=2)

# Save figure after cropping:
figname=prefix+str(frame)+'.eps'
plt.savefig(figname, format='eps', dpi=600)

plt.show()
