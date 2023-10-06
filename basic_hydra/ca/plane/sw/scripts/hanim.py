#!/usr/bin/env python

#=====================================================================
#   Animates h(x,0,t) for x > 0
#=====================================================================

#=====perform various generic imports=====
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.animation as anim
#=========================================

# Get the grid resolution:
file = open('src/parameters.f90').readlines()
for line in file:
   if 'ng=' in line:
      ng = int(line.split("=")[-1].strip())

N=ng*ng

# Open ecomp.asc file to get times:
in_file=open('evolution/ecomp.asc','r')
time, e1, e2, e3, e4, e5 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
time=np.array(time)

# Read in all data:
in_file=open('evolution/hh.r8','r')
raw_array=np.fromfile(in_file,dtype=np.float64)
in_file.close()

# For each frame, extract cross section:
nx=int(ng/2)
hc=np.empty([nx+1,len(time)])
for frame in range(len(time)):
   Z1=np.empty([ng,ng])
   Z1=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   hc[0:nx,frame]=Z1[nx:ng,nx]
   hc[nx,frame]=Z1[0,nx]

x=np.linspace(0.0,np.pi,nx+1)

# Create figure:
fig = plt.figure(figsize=[10.0,10.0])

# Add axes to figure
ax = fig.add_subplot(111,autoscale_on=False)

# Fix domain of view:
ax.set_xlim([0.0,np.pi])
ymax=abs(hc[0,0])/2.0
ax.set_ylim([-ymax,ymax])

# Formatting for displayed time:
time_template = 'time = %.3f'

# Set up list of images for animation:
ims=[]
# Repeat 1st frame 10 times:
for k in range(10):
   y=hc[:,0]

   # Create image of profile y = u(x,t):
   cims = ax.plot(x,y,c='k',lw=2)

   # Add time as a text string:
   cims.append(ax.text(0.02, 0.94, time_template%(0.0), transform=ax.transAxes, weight='bold', size=20))

   # Add to list of images containing all contours at each time:
   ims.append(cims)

# Add remaining frames:
for frame,t in enumerate(time):
   y=hc[:,frame]

   # Create image of profile y = u(x,t):
   cims = ax.plot(x,y,c='k',lw=2)

   # Add time as a text string:
   cims.append(ax.text(0.02, 0.94, time_template%(t), transform=ax.transAxes, weight='bold', size=20))

   # Add to list of images containing all contours at each time:
   ims.append(cims)

# Run animation:
ani = anim.ArtistAnimation(fig, ims, interval=100, repeat_delay=100, blit=False)

print
print(' Creating movie (evo.mp4) ... (this can take some time)')
ani.save('hc.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

print
print (' Display the movie by typing')
print
print (' mplayer hc.mp4 -loop 0')
print

plt.show()
