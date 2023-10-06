#! /usr/bin/python
import numpy as np
import matplotlib.pyplot as plt


# Open input file:
try:
   in_file=open('ene.dat','r')# try opening passed filename  
except IOError, message:# error if file not found 
   print >> sys.stderr, 'File could not be opened', message
   sys.exit()
  

#  Read in data and close input file:
raw_data = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

n=len(raw_data[0])
rocke=np.zeros(n-1)
rocpe=np.zeros(n-1)
ta=np.zeros(n-1)

t=raw_data[0]
ene=raw_data[1]
ke=raw_data[2]
pe=raw_data[3]
dkedt=raw_data[4]

#for i in range(0,n-1):
#  rocke[i]=(ke[i+1]-ke[i])/(t[i+1]-t[i])
#  rocpe[i]=-(pe[i+1]-pe[i])/(t[i+1]-t[i])
#  ta[i]=(t[i+1]+t[i])/2.0


fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(t,ke,'r-')
ax.plot(t,pe,'g-')
ax.plot(t,ene,'k-')

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.plot(t,dkedt,'k-')
#ax2.plot(t,pe,'g-')
#ax2.plot(ta,rocke,'r-')

plt.show()
