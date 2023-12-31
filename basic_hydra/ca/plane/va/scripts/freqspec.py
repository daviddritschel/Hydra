#======================================================================
# Creates frequency power spectra for the divergence sampled at a few
# specific points (nsamp x nsamp) from data generated by caps.f90.
#======================================================================

#=====perform various generic imports=====
import warnings,os,sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

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
#=========================================

# Get the sample size:
file = open('src/parameters.f90').readlines()
nsamp=''
for line in file:
   if 'nsamp=' in line:
      nsamp = line.split("=")[-1].strip()

n = int(nsamp)**2
nt = sum(1 for line in open('spectra/dsamp.asc'))
print(' There are ',nt,' times in the data.')

# Get hbar, cgw & cof:
for line in file:
   if 'hbar=' in line:
      hbar = float(line.split("=")[-1].strip().rstrip("d0"))
      cgw = float(line.split("=")[-2].strip().rstrip("d0,hbar"))
      cof = float(line.split("=")[-3].strip().rstrip("d0,cgw"))

# Compute buoyancy frequency N:
bvf=np.sqrt(3.0)*cgw/hbar
print(' N/f = ',bvf/cof)

twopi=2.0*np.pi
om1=cof/twopi
om2=bvf/twopi

time=np.zeros(nt)
d=np.zeros((nt,n))

in_file=open('spectra/dsamp.asc','r')
for j in range(nt):
   line=in_file.readline()  
   string=line.split()
   time[j]=float(string[0])
   for i in range(n):
      d[j,i]=float(string[i+1])

in_file.close()

print('')
print(' Creating frequency spectrum...')

# Create array holding the frequencies:
time_step = time[1] - time[0]
freqs = np.fft.fftfreq(time.size, time_step)

# Accumulate power spectra:
Sd=np.zeros(nt)
for k in range(n):
   ps=np.abs(np.fft.fft(d[:,k]))**2
   Sd+=ps
Sd=Sd/float(n)
   
# As signal is real, fold over and add spectrum:
nh_times=int((nt-1)/2)
for i in range(1,nh_times):
   Sd[i]=Sd[i]+Sd[nt-i]

logomega=np.log10(freqs[1:nh_times])
logSd=np.log10(Sd[1:nh_times])

#------------------------------------------------------------------------
# Plot S vs omega:
fig1 = plt.figure(1,figsize=[10,6])
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('$\log_{10}f$', fontsize=30)
ax1.set_ylabel('$\log_{10}{\mathcal{P}}_d$', fontsize=30)
ax1.axvline(np.log10(om1),color='r',linestyle='--')
ax1.axvline(np.log10(om2),color='b',linestyle='--')
# Limits in log_10{f}:
#xlims=[-2.5,1.5]
#ax1.set_xlim(xlims)
ax1.plot(logomega,logSd,c='k',lw=1)
fig1.savefig('d_fspec.eps', format='eps', dpi=600)

print()
print(' To view the frequency spectrum, type')
print()
print(' gv d_fspec.eps &')
print()
