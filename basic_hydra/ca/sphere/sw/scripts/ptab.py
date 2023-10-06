#!/usr/bin/env python

# This script plots data in table.txt in a chosen manner

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
mpl.rcParams['axes.linewidth'] = 1

#=================================================================
# Open data and read:
in_file=open('table.txt','r')
kd, rth, brms, tb, njet, angm = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
nd=len(kd)

#=================================================================
# Set up figure (generic parts):
fig = plt.figure(figsize=[8,8])  # Square figure
ax = plt.axes([0.2, 0.2, 0.7, 0.7])
#ax.set_aspect(1)

#=================================================================
# Select plot to make:
print()
print(' Choose one of the following plots to make:')
print()
print(' (1) N_J vs 1/tau for two selected values of b_rms')
print(' (2) N_J vs 1/tau for two selected values of t_b')
print(' (3) N_J vs b_rms for two selected values of tau')
print(' (4) N_J vs b_rms for two selected values of t_b')
print(' (5) N_J vs t_b for two selected values of b_rms')
print(' (6) N_J vs t_b for two selected values of tau')

print()
fopt_in = input(' Choice (default 1)? ')
fopt = int(fopt_in or 1)

if fopt==1:
   z = input(' 1st value of b_rms (default 0.01)? ')
   brms1 = float(z or 0.01)
   z = input(' 2nd value of b_rms (default 0.04)? ')
   brms2 = float(z or 0.04)
   z = input(' Value of t_b (default 1)? ')
   tb0 = float(z or 1)
   z = input(' Value of k_D (default 32)? ')
   kd0 = float(z or 32)

   ax.set_xlabel('$1/\\tau$', fontsize=30)
   ax.set_ylabel('$N_J$', fontsize=30)
   ax.set_title('$t_b = {x:.0f},$'.format(x=tb0)+'$\ \ \ \ k_{{\mathsf{{D}}}} = {x:.0f}$'.format(x=kd0), fontsize=30)
   ax.set_xscale('log')

   x1=[]
   y1=[]
   x2=[]
   y2=[]
   for k in range(nd):
      if abs(kd[k]-kd0) < 1.e-6:
         if abs(tb[k]-tb0) < 1.e-6:
            if abs(brms[k]-brms1) < 1.e-6:
               x1.append(rth[k])
               y1.append(njet[k])
            elif abs(brms[k]-brms2) < 1.e-6:
               x2.append(rth[k])
               y2.append(njet[k])

   n1=len(x1)
   if n1 > 0:
      x1=np.array(x1)
      y1=np.array(y1)
      # Sort so that x1 increases:
      ind=np.argsort(x1)
      x1s=np.empty(n1)
      y1s=np.empty(n1)
      for i in range(n1):
         x1s[i]=x1[ind[i]]
         y1s[i]=y1[ind[i]]
      ax.plot(x1s,y1s,c='k',lw=2,label='$b_{{\mathsf{{rms}}}} = {x:.2f}$'.format(x=brms1))
      ax.scatter(x1s, y1s, s=40, color='k', edgecolors='none')

   n2=len(x2)
   if n2 > 0:
      x2=np.array(x2)
      y2=np.array(y2)
      ind=np.argsort(x2)
      x2s=np.empty(n2)
      y2s=np.empty(n2)
      for i in range(n2):
         x2s[i]=x2[ind[i]]
         y2s[i]=y2[ind[i]]
      ax.plot(x2s,y2s,c='b',lw=2,label='$b_{{\mathsf{{rms}}}} = {x:.2f}$'.format(x=brms2))
      ax.scatter(x2s, y2s, s=40, color='b', edgecolors='none')

   ax.legend(loc='upper left',prop={'size':25})
   outfile='nj_vs_rth_for_tb'+str(tb0)+'kd'+str(kd0)+'.eps'

elif fopt==2:
   z = input(' 1st value of t_b (default 1)? ')
   tb1 = float(z or 1)
   z = input(' 2nd value of t_b (default 100)? ')
   tb2 = float(z or 100)
   z = input(' Value of b_rms (default 0.01)? ')
   brms0 = float(z or 0.01)
   z = input(' Value of k_D (default 32)? ')
   kd0 = float(z or 32)

   ax.set_xlabel('$1/\\tau$', fontsize=30)
   ax.set_ylabel('$N_J$', fontsize=30)
   ax.set_title('$b_{{\mathsf{{rms}}}} = {x:.2f},$'.format(x=brms0)+'$\ \ \ \ k_{{\mathsf{{D}}}} = {x:.0f}$'.format(x=kd0), fontsize=30)
   ax.set_xscale('log')

   x1=[]
   y1=[]
   x2=[]
   y2=[]
   for k in range(nd):
      if abs(kd[k]-kd0) < 1.e-6:
         if abs(brms[k]-brms0) < 1.e-6:
            if abs(tb[k]-tb1) < 1.e-6:
               x1.append(rth[k])
               y1.append(njet[k])
            elif abs(tb[k]-tb2) < 1.e-6:
               x2.append(rth[k])
               y2.append(njet[k])


   n1=len(x1)
   if n1 > 0:
      x1=np.array(x1)
      y1=np.array(y1)
      # Sort so that x1 increases:
      ind=np.argsort(x1)
      x1s=np.empty(n1)
      y1s=np.empty(n1)
      for i in range(n1):
         x1s[i]=x1[ind[i]]
         y1s[i]=y1[ind[i]]
      ax.plot(x1s,y1s,c='k',lw=2,label='$t_b = {x:.0f}$'.format(x=tb1))
      ax.scatter(x1s, y1s, s=40, color='k', edgecolors='none')

   n2=len(x2)
   if n2 > 0:
      x2=np.array(x2)
      y2=np.array(y2)
      ind=np.argsort(x2)
      x2s=np.empty(n2)
      y2s=np.empty(n2)
      for i in range(n2):
         x2s[i]=x2[ind[i]]
         y2s[i]=y2[ind[i]]
      ax.plot(x2s,y2s,c='b',lw=2,label='$t_b = {x:.0f}$'.format(x=tb2))
      ax.scatter(x2s, y2s, s=40, color='b', edgecolors='none')
   
   ax.legend(loc='upper left',prop={'size':25})
   outfile='nj_vs_rth_for_brms'+str(brms0)+'kd'+str(kd0)+'.eps'

elif fopt==3:
   z = input(' 1st value of tau (default 100)? ')
   tau1 = float(z or 100)
   rth1 = 1.0/tau1
   z = input(' 2nd value of tau (default 1000)? ')
   tau2 = float(z or 1000)
   rth2 = 1.0/tau2
   z = input(' Value of t_b (default 1)? ')
   tb0 = float(z or 1)
   z = input(' Value of k_D (default 32)? ')
   kd0 = float(z or 32)

   ax.set_xlabel('$b_{\mathsf{rms}}$', fontsize=30)
   ax.set_ylabel('$N_J$', fontsize=30)
   ax.set_title('$t_b = {x:.0f},$'.format(x=tb0)+'$\ \ \ \ k_{{\mathsf{{D}}}} = {x:.0f}$'.format(x=kd0), fontsize=30)

   x1=[]
   y1=[]
   x2=[]
   y2=[]
   for k in range(nd):
      if abs(kd[k]-kd0) < 1.e-6:
         if abs(tb[k]-tb0) < 1.e-6:
            if abs(rth[k]-rth1) < 1.e-6:
               x1.append(brms[k])
               y1.append(njet[k])
            elif abs(rth[k]-rth2) < 1.e-6:
               x2.append(brms[k])
               y2.append(njet[k])

   n1=len(x1)
   if n1 > 0:
      x1=np.array(x1)
      y1=np.array(y1)
      # Sort so that x1 increases:
      ind=np.argsort(x1)
      x1s=np.empty(n1)
      y1s=np.empty(n1)
      for i in range(n1):
         x1s[i]=x1[ind[i]]
         y1s[i]=y1[ind[i]]
      ax.plot(x1s,y1s,c='k',lw=2,label='$\\tau = {x:.0f}$'.format(x=tau1))
      ax.scatter(x1s, y1s, s=40, color='k', edgecolors='none')

   n2=len(x2)
   if n2 > 0:
      x2=np.array(x2)
      y2=np.array(y2)
      ind=np.argsort(x2)
      x2s=np.empty(n2)
      y2s=np.empty(n2)
      for i in range(n2):
         x2s[i]=x2[ind[i]]
         y2s[i]=y2[ind[i]]
      ax.plot(x2s,y2s,c='b',lw=2,label='$\\tau = {x:.0f}$'.format(x=tau2))
      ax.scatter(x2s, y2s, s=40, color='b', edgecolors='none')

   ax.legend(loc='upper right',prop={'size':25})
   outfile='nj_vs_brms_for_tb'+str(tb0)+'kd'+str(kd0)+'.eps'

elif fopt==4:
   z = input(' 1st value of t_b (default 1)? ')
   tb1 = float(z or 1)
   z = input(' 2nd value of t_b (default 100)? ')
   tb2 = float(z or 100)
   z = input(' Value of tau (default 1000)? ')
   tau0 = float(z or 1000)
   rth0 = 1.0/tau0
   z = input(' Value of k_D (default 32)? ')
   kd0 = float(z or 32)

   ax.set_xlabel('$b_{\mathsf{rms}}$', fontsize=30)
   ax.set_ylabel('$N_J$', fontsize=30)
   ax.set_title('$\\tau = {x:.0f},$'.format(x=tau0)+'$\ \ \ \ k_{{\mathsf{{D}}}} = {x:.0f}$'.format(x=kd0), fontsize=30)

   x1=[]
   y1=[]
   x2=[]
   y2=[]
   for k in range(nd):
      if abs(kd[k]-kd0) < 1.e-6:
         if abs(rth[k]-rth0) < 1.e-6:
            if abs(tb[k]-tb1) < 1.e-6:
               x1.append(brms[k])
               y1.append(njet[k])
            elif abs(tb[k]-tb2) < 1.e-6:
               x2.append(brms[k])
               y2.append(njet[k])

   n1=len(x1)
   if n1 > 0:
      x1=np.array(x1)
      y1=np.array(y1)
      # Sort so that x1 increases:
      ind=np.argsort(x1)
      x1s=np.empty(n1)
      y1s=np.empty(n1)
      for i in range(n1):
         x1s[i]=x1[ind[i]]
         y1s[i]=y1[ind[i]]
      ax.plot(x1s,y1s,c='k',lw=2,label='$t_b = {x:.0f}$'.format(x=tb1))
      ax.scatter(x1s, y1s, s=40, color='k', edgecolors='none')

   n2=len(x2)
   if n2 > 0:
      x2=np.array(x2)
      y2=np.array(y2)
      ind=np.argsort(x2)
      x2s=np.empty(n2)
      y2s=np.empty(n2)
      for i in range(n2):
         x2s[i]=x2[ind[i]]
         y2s[i]=y2[ind[i]]
      ax.plot(x2s,y2s,c='b',lw=2,label='$t_b = {x:.0f}$'.format(x=tb2))
      ax.scatter(x2s, y2s, s=40, color='b', edgecolors='none')

   ax.legend(loc='upper right',prop={'size':25})
   outfile='nj_vs_brms_for_tau'+str(tau0)+'kd'+str(kd0)+'.eps'

elif fopt==5:
   z = input(' 1st value of b_rms (default 0.01)? ')
   brms1 = float(z or 0.01)
   z = input(' 2nd value of b_rms (default 0.04)? ')
   brms2 = float(z or 0.04)
   z = input(' Value of tau (default 1000)? ')
   tau0 = float(z or 1000)
   rth0 = 1.0/tau0
   z = input(' Value of k_D (default 32)? ')
   kd0 = float(z or 32)

   ax.set_xlabel('$t_b$', fontsize=30)
   ax.set_ylabel('$N_J$', fontsize=30)
   ax.set_title('$\\tau = {x:.0f},$'.format(x=tau0)+'$\ \ \ \ k_{{\mathsf{{D}}}} = {x:.0f}$'.format(x=kd0), fontsize=30)
   ax.set_xscale('log')

   x1=[]
   y1=[]
   x2=[]
   y2=[]
   for k in range(nd):
      if abs(kd[k]-kd0) < 1.e-6:
         if abs(rth[k]-rth0) < 1.e-6:
            if abs(brms[k]-brms1) < 1.e-6:
               x1.append(tb[k])
               y1.append(njet[k])
            elif abs(brms[k]-brms2) < 1.e-6:
               x2.append(tb[k])
               y2.append(njet[k])

   n1=len(x1)
   if n1 > 0:
      x1=np.array(x1)
      y1=np.array(y1)
      # Sort so that x1 increases:
      ind=np.argsort(x1)
      x1s=np.empty(n1)
      y1s=np.empty(n1)
      for i in range(n1):
         x1s[i]=x1[ind[i]]
         y1s[i]=y1[ind[i]]
      ax.plot(x1s,y1s,c='k',lw=2,label='$b_{{\mathsf{{rms}}}} = {x:.2f}$'.format(x=brms1))
      ax.scatter(x1s, y1s, s=40, color='k', edgecolors='none')

   n2=len(x2)
   if n2 > 0:
      x2=np.array(x2)
      y2=np.array(y2)
      ind=np.argsort(x2)
      x2s=np.empty(n2)
      y2s=np.empty(n2)
      for i in range(n2):
         x2s[i]=x2[ind[i]]
         y2s[i]=y2[ind[i]]
      ax.plot(x2s,y2s,c='b',lw=2,label='$b_{{\mathsf{{rms}}}} = {x:.2f}$'.format(x=brms2))
      ax.scatter(x2s, y2s, s=40, color='b', edgecolors='none')

   ax.legend(loc='center right',prop={'size':25})
   outfile='nj_vs_tb_for_tau'+str(tau0)+'kd'+str(kd0)+'.eps'

elif fopt==6:
   z = input(' 1st value of tau (default 100)? ')
   tau1 = float(z or 100)
   rth1 = 1.0/tau1
   z = input(' 2nd value of tau (default 1000)? ')
   tau2 = float(z or 1000)
   rth2 = 1.0/tau2
   z = input(' Value of b_rms (default 0.01)? ')
   brms0 = float(z or 0.01)
   z = input(' Value of k_D (default 32)? ')
   kd0 = float(z or 32)

   ax.set_xlabel('$t_b$', fontsize=30)
   ax.set_ylabel('$N_J$', fontsize=30)
   ax.set_title('$b_{{\mathsf{{rms}}}} = {x:.2f},$'.format(x=brms0)+'$\ \ \ \ k_{{\mathsf{{D}}}} = {x:.0f}$'.format(x=kd0), fontsize=30)
   ax.set_xscale('log')

   x1=[]
   y1=[]
   x2=[]
   y2=[]
   for k in range(nd):
      if abs(kd[k]-kd0) < 1.e-6:
         if abs(brms[k]-brms0) < 1.e-6:
            if abs(rth[k]-rth1) < 1.e-6:
               x1.append(tb[k])
               y1.append(njet[k])
            elif abs(rth[k]-rth2) < 1.e-6:
               x2.append(tb[k])
               y2.append(njet[k])

   n1=len(x1)
   if n1 > 0:
      x1=np.array(x1)
      y1=np.array(y1)
      # Sort so that x1 increases:
      ind=np.argsort(x1)
      x1s=np.empty(n1)
      y1s=np.empty(n1)
      for i in range(n1):
         x1s[i]=x1[ind[i]]
         y1s[i]=y1[ind[i]]
      ax.plot(x1s,y1s,c='k',lw=2,label='$\\tau = {x:.0f}$'.format(x=tau1))
      ax.scatter(x1s, y1s, s=40, color='k', edgecolors='none')

   n2=len(x2)
   if n2 > 0:
      x2=np.array(x2)
      y2=np.array(y2)
      ind=np.argsort(x2)
      x2s=np.empty(n2)
      y2s=np.empty(n2)
      for i in range(n2):
         x2s[i]=x2[ind[i]]
         y2s[i]=y2[ind[i]]
      ax.plot(x2s,y2s,c='b',lw=2,label='$\\tau = {x:.0f}$'.format(x=tau2))
      ax.scatter(x2s, y2s, s=40, color='b', edgecolors='none')

   ax.legend(loc='lower right',prop={'size':25})
   outfile='nj_vs_tb_for_brms'+str(brms0)+'kd'+str(kd0)+'.eps'

#=========================================================================

fig.savefig(outfile, format='eps', dpi=300)
print()
print(' To view the image, type')
print()
print(' gv ',outfile,' &')
print()
