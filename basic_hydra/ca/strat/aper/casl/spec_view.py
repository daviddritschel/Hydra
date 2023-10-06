#! /usr/bin/python
import subprocess as sbpc
import sys
sys.path.append('/user/stuart/scripts/modules')
import numpy as np
import matplotlib.pyplot as plt
import argparse


def parse_args():
  parser = argparse.ArgumentParser(prog='framepy.py')
  parser.add_argument('input', metavar='input_filename' , type=str , help='Input file')
  parser.add_argument('-xlims', metavar='x_lim' , type=float , nargs=2, default='0.0 3.0'.split(),help='Lower and upper limits for x range')
  parser.add_argument('-ylims', metavar='y_lim' , type=float , nargs=2, default='-13.0 3.0'.split(),help='Lower and upper limits for y range')
  args = parser.parse_args() 
  return args

#Define a routine for adding zeros to the frame output files
def add_nulls(num, cnt=3):
  cnt = cnt - len(str(num))
  nulls = '0' * cnt
  return '%s%s' % (nulls, num)

def running(args):

  # Open input file:
  try:
     in_file=open(args.input,'r')# try opening passed filename  
  except IOError, message:# error if file not found 
     print >> sys.stderr, 'File could not be opened', message
     sys.exit()


  
#  Read in data and close input file:
  raw_data = np.array(np.fromfile(file=in_file,dtype=float,sep='\n'))
  in_file.close()

#  t0=raw_data[0]   
#  sumzpsec=raw_data[1]   
  kmax=int(raw_data[3])   

# Set the number of frames:
  nframes = int(len(raw_data)/(2*kmax+3))  
  print nframes
  print kmax

# Shape the data array into a useful shape for plotting:
  frames=range(0,nframes)
  tim_eles = [i*(2*kmax+4)+j for i in frames for j in range(0,4)]
#  print tim_eles
  shaped_data = np.delete(raw_data,tim_eles)[0:(2*kmax+4)*nframes].reshape((nframes,kmax,2))
  global ic 
  # Grab the correct sub-array for plotting  
  ic = 0
  
  fig = plt.figure(1)
  ax = fig.add_subplot(111)
  im = ax.plot(shaped_data[ic].transpose()[0],shaped_data[ic].transpose()[1],'k-')
  ax.set_title('Frame no.: '+str(ic+1))
  xlimits=[float(x) for x in args.xlims]
  ylimits=[float(y) for y in args.ylims]
  ax.set_xlim(xlimits)  
  ax.set_ylim(ylimits)

  def on_press(event):
    global ic
    if event.key=='=':
      axes=event.canvas.figure.get_axes()[0]
      if ic != nframes-1:
        ic+=1
      xlim=axes.get_xlim()
      ylim=axes.get_ylim()
      axes.clear()
      axes.plot(shaped_data[ic].transpose()[0],shaped_data[ic].transpose()[1],'k-')
      axes.set_title('Frame no.: '+str(ic+1))
      axes.set_xlim(xlim)  
      axes.set_ylim(ylim)
      plt.draw()
    if event.key=='-':
      axes=event.canvas.figure.get_axes()[0]
      if ic != 0:
        ic-=1
      xlim=axes.get_xlim()
      ylim=axes.get_ylim()
      axes.clear()
      axes.plot(shaped_data[ic].transpose()[0],shaped_data[ic].transpose()[1],'k-')
      axes.set_title('Frame no.: '+str(ic+1))
      axes.set_xlim(xlim)  
      axes.set_ylim(ylim)
      plt.draw()


#    print axes[0].get_title()
#    print 'you pressed', event.button, event.xdata, event.ydata
  cid=fig.canvas.mpl_connect('key_press_event',on_press)
  plt.show()
 

if __name__ == '__main__':
  args = parse_args()
  running(args)



