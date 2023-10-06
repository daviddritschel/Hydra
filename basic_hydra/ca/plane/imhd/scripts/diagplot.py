#=================================================================
#    Plots diagnostics in ene.asc, norms.asc and circ.asc

#       Written 4/9/2017 by D G Dritschel @ St Andrews
#=================================================================

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

## global settings

# set tick label size:
label_size = 20
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

#------------------------------------------------------------
#Get delta and gamma from job_info file for scaling purposes:
f = open('job_info','r')
for j in range(28):
   line=f.readline()
string=f.readline().split()
delta='%.5f'%np.float(string[-1])
string=f.readline().split()
gamma='%.3f'%np.float(string[-1])
f.close()

print(' delta = ',delta,'  gamma = ',gamma)

#------------------------------------------------------------
#Create figures:

#Line colours:
col=['k','b','r','m','g','c']

title_in=input('Add delta and gamma to title? (default y) ')
title=str(title_in or 'y')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#energies vs t
fig1 = plt.figure(1,figsize=[8,8])
ax1 = plt.axes([0.16, 0.16, 0.72, 0.72])
ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('Energy', fontsize=30)
if title=='y':
   ax1.set_title('$\\gamma=$'+str(gamma)+'    $\\delta=$'+str(delta), fontsize=30)

fname='ene.asc'
n_lines = sum(1 for line in open(fname))
time=np.zeros(n_lines)
eneu=np.zeros(n_lines)
eneb=np.zeros(n_lines)
enet=np.zeros(n_lines)
diag_file = open(fname,'r')
for j in range(n_lines):
   string=diag_file.readline().split()
   time[j]=float(string[0])
   eneu[j]=float(string[1])
   eneb[j]=float(string[2])
   enet[j]=float(string[3])
diag_file.close

ax1.plot(time, eneu, c=col[0], lw=2, label='$E_{\mathsf{u}}$')
ax1.plot(time, eneb, c=col[1], lw=2, label='$E_{\mathsf{b}}$')
ax1.plot(time, enet, c=col[2], lw=2, label='$E_{\mathsf{u}}+E_{\mathsf{b}}$')

ax1.legend(loc='upper right', fontsize=20)

fig1.savefig('ene.png', bbox_inches='tight', pad_inches = 0.025, dpi=100)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#omega^2 and j^2 vs t
fig2 = plt.figure(2,figsize=[8,8])
ax2 = plt.axes([0.16, 0.16, 0.72, 0.72])
ax2.set_xlabel('$t$', fontsize=30)
ax2.set_ylabel('Enstrophy', fontsize=30)
if title=='y':
   ax2.set_title('$\\gamma=$'+str(gamma)+'    $\\delta=$'+str(delta), fontsize=30)

fname='norms.asc'
n_lines = sum(1 for line in open(fname))
time=np.zeros(n_lines)
ensu=np.zeros(n_lines)
ensb=np.zeros(n_lines)
diag_file = open(fname,'r')
for j in range(n_lines):
   string=diag_file.readline().split()
   time[j]=float(string[0])
   ensu[j]=float(string[1])
   ensb[j]=float(string[2])
diag_file.close

ax2.plot(time, ensu, c=col[0], lw=2, label='$\\langle{\\omega^2}\\rangle$')
ax2.plot(time, ensb, c=col[1], lw=2, label='$\\langle{j^2}\\rangle$')

ax2.legend(loc='upper left', fontsize=20)

fig2.savefig('norms.png', bbox_inches='tight', pad_inches = 0.025, dpi=100)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#circulation (Gamma) vs t  and  dGamma/dt vs t
fig3 = plt.figure(3,figsize=[8,8])
ax3 = plt.axes([0.16, 0.16, 0.72, 0.72])
ax3.set_xlabel('$t$', fontsize=30)
ax3.set_ylabel('$\\Gamma$', fontsize=30)
if title=='y':
   ax3.set_title('$\\gamma=$'+str(gamma)+'    $\\delta=$'+str(delta), fontsize=30)

fig4 = plt.figure(4,figsize=[8,8])
ax4 = plt.axes([0.16, 0.16, 0.72, 0.72])
ax4.set_xlabel('$t$', fontsize=30)
ax4.set_ylabel('$d\\Gamma/dt$', fontsize=30)
if title=='y':
   ax4.set_title('$\\gamma=$'+str(gamma)+'    $\\delta=$'+str(delta), fontsize=30)

fname='circulation.asc'
n_lines = sum(1 for line in open(fname))
time=np.zeros(n_lines)
circ=np.zeros(n_lines)
dcdt=np.zeros(n_lines)
diag_file = open(fname,'r')
for j in range(n_lines):
   string=diag_file.readline().split()
   time[j]=float(string[0])
   circ[j]=float(string[1])
   dcdt[j]=float(string[3])
diag_file.close

ax3.plot(time, circ, c=col[0], lw=2)
ax4.plot(time, dcdt, c=col[0], lw=2)

fig3.savefig('circ.png', bbox_inches='tight', pad_inches = 0.025, dpi=100)
fig4.savefig('dcdt.png', bbox_inches='tight', pad_inches = 0.025, dpi=100)

plt.show()
