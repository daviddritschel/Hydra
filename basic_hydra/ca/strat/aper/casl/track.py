#! /usr/bin/python
import numpy as np
import matplotlib.pyplot as plt


# Open input file:
try:
   in_file=open('bb.dat','r')# try opening passed filename  
except IOError, message:# error if file not found 
   print >> sys.stderr, 'File could not be opened', message
   sys.exit()
  

#  Set up grid size etc.:
nx = 1024
ny = 128

ellx = 8.0
hlx = ellx/2.0
delx = ellx/float(nx)

bmid=-0.5

#  Read in data and close input file:
raw_data = np.array(np.fromfile(file=in_file,dtype=float,sep='\n'))
in_file.close()

# Set the number of frames:
nframes = int(len(raw_data)/((nx+1)*(ny+1)+1))  


# Shape the data array into a useful shape for plotting:
frames=range(0,nframes)
tim_eles = [i*((nx+1)*(ny+1)+1) for i in frames]
shaped_data = np.delete(raw_data,tim_eles)[0:((nx+1)*(ny+1)*nframes)].reshape((nframes,nx+1,ny+1))

x=[]
t=[]
x2=[]
for frame in frames:
  t.append(raw_data[tim_eles[frame]])
  Z=shaped_data[frame]
  for ix in range(0,nx+1):
    if Z[ix,0] > bmid:
      p=(Z[ix-1,0]-bmid)/(Z[ix-1,0]-Z[ix,0])
      x.append((ix-1+p)*delx-hlx)
      break    
  for ix in range(nx,0,-1):
    if Z[ix,ny] < bmid:
      p=(Z[ix,ny]-bmid)/(Z[ix,ny]-Z[ix+1,ny])
      x2.append(hlx-(ix+p)*delx)
      break    


a=np.arange(0,10,0.01)
b=0.4*a
b2=0.5*a
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(t,x,'r+-')
ax.plot(t,x2,'bx-')
ax.plot(a,b,'k-')
ax.plot(a,b2,'b-')

#im = ax.imshow(Z.transpose(),cmap=plt.cm.jet,origin='lower',interpolation='nearest')

plt.show()
