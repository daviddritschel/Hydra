#!/usr/bin/env python

# This script plots a chosen field (either full, balanced or imbalanced)
# for data in 4 separate directories, as indicated below.  Make sure to
# adjust the titles for each directory.

#========== Perform the generic imports =========
import warnings
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

#====================== Function definitions =======================
def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

    emag=1.0
    rmult=fmax-fmin
    while rmult < 10:
       emag=emag/10
       rmult=rmult*10

    while rmult >= 100:
       emag=emag*10
       rmult=rmult/10

    kmult=int(rmult/10)

    if kmult < 1:
       ci=emag
    elif kmult < 2:
       ci=2*emag
    elif kmult < 4:
       ci=4*emag
    elif kmult < 8:
       ci=10*emag
    else:
       ci=20*emag

    return ci

#=================================================================
# Specify the data directories (need final /):
basedir='/home/dgd/data/hydra/ca/plane/'
dir_list=[basedir+'sw/caps/ranpv/ng256kd6dr10r001/', \
          basedir+'va/caps/ranpv/ng256H0.1kd6dr10r001/', \
          basedir+'va/caps/ranpv/ng256H0.2kd6dr10r001/', \
          basedir+'va/caps/ranpv/ng256H0.4kd6dr10r001/']
# Corresponding titles (supplemented by chosen field):
label_list=['$H=0$','$H=0.1$','$H=0.2$','$H=0.4$']

# Specify grid resolution (ng x ng):
ng=256
N=ng*ng

#=================================================================
# Select data to compare:
print()
print(' Which field do you wish to show:')
print()
print(' (1) h_tilde')
print(' (2) zeta')
print(' (3) delta')
print(' (4) gamma')
print(' (5) gamma_tilde')
print(' (6) q or')
print(' (7) P')

print()
k_in = input(' Option (default 3)? ')
k = int(k_in or 3)

if k==1:
   dfile='hh.r8'
   field='h'
elif k==2:
   dfile='zz.r8'
   field='z'
elif k==3:
   dfile='dd.r8'
   field='d'
elif k==4:
   dfile='gg.r8'
   field='g'
elif k==5:
   dfile='gt.r8'
   field='gt'
elif k==6:
   dfile='qq.r8'
   field='q'
elif k==7:
   dfile='pp.r8'
   field='p'
else:
   print(' *** Not a valid option!  Stopping! ***')
   exit

#=================================================================
# Select full, balanced or imbalanced fields:
opt_list=[]
stitle_list=[]
for dir in dir_list:
   print(' For data in',dir)
   opt_in = input(' show (1) full, (2) balanced or (3) imbalanced field (default 1)? ')
   option = int(opt_in or 1)
   if k==1:
      if option==1:
         stitle='$\\tilde{h}$'
      elif option==2:
         stitle='$\\tilde{h}_{\mathsf{b}}$'
      else:
         stitle='$\\tilde{h}_{\mathsf{i}}$'
   elif k==2:
      if option==1:
         stitle='$\\zeta$'
      elif option==2:
         stitle='$\\zeta_{\mathsf{b}}$'
      else:
         stitle='$\\zeta_{\mathsf{i}}$'
   elif k==3:
      if option==1:
         stitle='$\\delta$'
      elif option==2:
         stitle='$\\delta_{\mathsf{b}}$'
      else:
         stitle='$\\delta_{\mathsf{i}}$'
   elif k==4:
      if option==1:
         stitle='$\\gamma$'
      elif option==2:
         stitle='$\\gamma_{\mathsf{b}}$'
      else:
         stitle='$\\gamma_{\mathsf{i}}$'
   elif k==5:
      if option==1:
         stitle='$\\tilde\\gamma$'
      elif option==2:
         stitle='$\\tilde\\gamma_{\mathsf{b}}$'
      else:
         stitle='$\\tilde\\gamma_{\mathsf{i}}$'
   elif k==6:
      if option==1:
         stitle='$q$'
      elif option==2:
         stitle='$q_{\mathsf{b}}$'
      else:
         stitle='$q_{\mathsf{i}}$'
   else:
      if option==1:
         stitle='$P$'
      elif option==2:
         stitle='$P_{\mathsf{b}}$'
      else:
         stitle='$P_{\mathsf{i}}$'

   opt_list.append(option-1)
   stitle_list.append(stitle)
   print()

# Label suffixes to indicate full, balanced or imbalanced fields:
label_suffix=['','(b)','(i)']
# Prefix to add to datafile:
prefile=['','b','']

# Select time to show:
t_in = input(' Time to show (default 500)? ')
t = float(t_in or 500.0)

# Output filename:
outfile=field+'_n'+str(ng)+'_t'+str(int(t+0.01))+'.eps'

#=================================================================
# Open ecomp.asc file in one directory to get time between frames:
in_file=open(dir_list[0]+'evolution/ecomp.asc','r')
time, e1, e2, e3, e4, e5 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

dt=time[1]-time[0]
# Frame corresponding to time chosen:
frame=int((t+0.0001)/dt)

#==============================================================================
# Set up figure:
fig, (ax1, ax2, ax3, ax4) = plt.subplots(figsize=[20,5], nrows=1, ncols=4)

# Customise tick values:
ax1.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax4.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax4.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

ax1.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax4.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax4.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

#=================================================================
# Read data into arrays for plotting:
k=0
Z1=np.empty([ng+1,ng+1])
in_file=open(dir_list[k]+'evolution/'+prefile[opt_list[k]]+dfile,'r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

if opt_list[k] == 2:
   Z=np.empty([ng+1,ng+1])
   in_file=open(dir_list[k]+'evolution/b'+dfile,'r')
   raw_array=np.fromfile(in_file,dtype=np.float64)
   in_file.close()
   Z[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
   Z1=Z1-Z

k=1
Z2=np.empty([ng+1,ng+1])
in_file=open(dir_list[k]+'evolution/'+prefile[opt_list[k]]+dfile,'r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

if opt_list[k] == 2:
   Z=np.empty([ng+1,ng+1])
   in_file=open(dir_list[k]+'evolution/b'+dfile,'r')
   raw_array=np.fromfile(in_file,dtype=np.float64)
   in_file.close()
   Z[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
   Z2=Z2-Z

k=2
Z3=np.empty([ng+1,ng+1])
in_file=open(dir_list[k]+'evolution/'+prefile[opt_list[k]]+dfile,'r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

if opt_list[k] == 2:
   Z=np.empty([ng+1,ng+1])
   in_file=open(dir_list[k]+'evolution/b'+dfile,'r')
   raw_array=np.fromfile(in_file,dtype=np.float64)
   in_file.close()
   Z[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
   Z3=Z3-Z

k=3
Z4=np.empty([ng+1,ng+1])
in_file=open(dir_list[k]+'evolution/'+prefile[opt_list[k]]+dfile,'r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()
Z4[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

if opt_list[k] == 2:
   Z=np.empty([ng+1,ng+1])
   in_file=open(dir_list[k]+'evolution/b'+dfile,'r')
   raw_array=np.fromfile(in_file,dtype=np.float64)
   in_file.close()
   Z[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
   Z4=Z4-Z

# Add periodic edges:
Z1[ng,0:ng]=Z1[0,0:ng]
Z1[0:ng+1,ng]=Z1[0:ng+1,0]
Z2[ng,0:ng]=Z2[0,0:ng]
Z2[0:ng+1,ng]=Z2[0:ng+1,0]
Z3[ng,0:ng]=Z3[0,0:ng]
Z3[0:ng+1,ng]=Z3[0:ng+1,0]
Z4[ng,0:ng]=Z4[0,0:ng]
Z4[0:ng+1,ng]=Z4[0:ng+1,0]

# Work out the overall min/max values:
zmin1=np.amin(Z1)
zmin2=np.amin(Z2)
zmin3=np.amin(Z3)
zmin4=np.amin(Z4)
zmax1=np.amax(Z1)
zmax2=np.amax(Z2)
zmax3=np.amax(Z3)
zmax4=np.amax(Z4)

print()
print(' Minimum and maximum field values for data in each directory:')
print()
print('  Directory    Min field value      Max field value')
print(' -----------   ---------------      ---------------')
print('     1st       ',"{:14.10f}".format(zmin1),'     ',"{:14.10f}".format(zmax1))
print('     2nd       ',"{:14.10f}".format(zmin2),'     ',"{:14.10f}".format(zmax2))
print('     3rd       ',"{:14.10f}".format(zmin3),'     ',"{:14.10f}".format(zmax3))
print('     4th       ',"{:14.10f}".format(zmin4),'     ',"{:14.10f}".format(zmax4))

# Use different limits for each simulation:
zmag=max(abs(zmin1),zmax1)
zmin1=-zmag
zmax1= zmag

zmag=max(abs(zmin2),zmax2)
zmin2=-zmag
zmax2= zmag

zmag=max(abs(zmin3),zmax3)
zmin3=-zmag
zmax3= zmag

zmag=max(abs(zmin4),zmax4)
zmin4=-zmag
zmax4= zmag

# Obtain contour levels for plotting the colorbars:
dz=contint(zmin1,zmax1)
jmin=-int(-zmin1/dz)
jmax= int( zmax1/dz)
clevels1=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin2,zmax2)
jmin=-int(-zmin2/dz)
jmax= int( zmax2/dz)
clevels2=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin3,zmax3)
jmin=-int(-zmin3/dz)
jmax= int( zmax3/dz)
clevels3=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin4,zmax4)
jmin=-int(-zmin4/dz)
jmax= int( zmax4/dz)
clevels4=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

# Plot the images in an array with individual colourbars:
im1=ax1.imshow(Z1,cmap=cm.seismic,vmin=zmin1,vmax=zmax1,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im1, cax=cax, ticks=clevels1)
#setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=16)

im2=ax2.imshow(Z2,cmap=cm.seismic,vmin=zmin2,vmax=zmax2,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im2, cax=cax, ticks=clevels2)
#setp(cbar.ax.yaxis.set_ticklabels(clevels2), fontsize=16)

im3=ax3.imshow(Z3,cmap=cm.seismic,vmin=zmin3,vmax=zmax3,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im3, cax=cax, ticks=clevels3)
#setp(cbar.ax.yaxis.set_ticklabels(clevels3), fontsize=16)

im4=ax4.imshow(Z4,cmap=cm.seismic,vmin=zmin4,vmax=zmax4,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax4)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im4, cax=cax, ticks=clevels4)
#setp(cbar.ax.yaxis.set_ticklabels(clevels4), fontsize=16)

# Add labels for each simulation:
ax1.set_title(stitle_list[0]+',  '+label_list[0], fontsize=20)
ax2.set_title(stitle_list[1]+',  '+label_list[1], fontsize=20)
ax3.set_title(stitle_list[2]+',  '+label_list[2], fontsize=20)
ax4.set_title(stitle_list[3]+',  '+label_list[3], fontsize=20)

ax1.label_outer()
ax2.label_outer()
ax3.label_outer()
ax4.label_outer()

# Fine-tune figure; hide y ticks for right plots
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)

#fig.subplots_adjust(wspace=0.8, hspace=0.5, top=0.6, bottom=0.05)
fig.subplots_adjust(wspace=0.8, hspace=-0.1)

#=========================================================================
# Save image:
fig.savefig(outfile, format='eps', dpi=300)

print()
print(' To view the image, type')
print()
print(' gv ',outfile)
print()
