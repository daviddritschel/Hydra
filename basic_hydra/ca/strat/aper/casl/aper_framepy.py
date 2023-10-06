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
  parser.add_argument('-frame', metavar='frame_no' , type=int,default='1', help='Frame number to examine')
  parser.add_argument('-ndim', metavar='n_dim' , type=int , nargs=2, default='1024 128'.split(), help='Number of x, and y grid points: nx ny')
  parser.add_argument('-lims', metavar='v_lim' , type=float , nargs=2, help='Lower and upper limits for the colourmap range')
  parser.add_argument('-xlims', metavar='xlims' , type=int , nargs=2, help='Lower and upper x grid points')
  parser.add_argument('-ylims', metavar='ylims' , type=int , nargs=2, help='Lower and upper y grid points')
  parser.add_argument('-cb', action='store_true' , help='Add colourbar to plots')  
  parser.add_argument('-noticks', action='store_true' , help='Remove tickmark numbers from the axes')  
  parser.add_argument('-notitle', action='store_true' , help='Remove title with frame number')  
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
  

#  Set up grid point numbers
  nx = int(args.ndim[0])
  ny = int(args.ndim[1])

#  Read in data and close input file:
  raw_data = np.array(np.fromfile(file=in_file,dtype=float,sep='\n'))
  in_file.close()

# Set the number of frames:
  nframes = int(len(raw_data)/((nx+1)*(ny+1)+1))  
  print nframes
  if args.frame > nframes:
    print >> sys.stderr, 'Cannot find that many frames' 
    sys.exit()

# Shape the data array into a useful shape for plotting:
  frames=range(0,nframes)
  tim_eles = [i*((nx+1)*(ny+1)+1) for i in frames]
  shaped_data = np.delete(raw_data,tim_eles)[0:((nx+1)*(ny+1)*nframes)].reshape((nframes,nx+1,ny+1))
  global ic 
  # Grab the correct sub-array for plotting  
  ic = args.frame-1
  Z=shaped_data[(args.frame-1)]

  print'Frame '+str(args.frame)+' min/max: %8.5f , %8.5f' % (Z.min(),Z.max())
  print 'Overall data min/max: %8.5f , %8.5f' % (shaped_data.min(),shaped_data.max())

#  lev_min=shaped_data.min()
#  lev_max=shaped_data.max()    
  lev_min=Z.min()
  lev_max=Z.max()    

  if args.lims:
    lev_min=args.lims[0]
    lev_max=args.lims[1]


  fig = plt.figure(1)
  ax = fig.add_subplot(111)
	
  im = ax.imshow(Z.transpose(),cmap=plt.cm.jet,vmin=lev_min,vmax=lev_max,origin='lower',interpolation='nearest')

  if args.xlims:
    ax.set_xlim(args.xlims)  
  if args.ylims:
    ax.set_ylim(args.ylims)


# Set ticks or not:
  if args.noticks:
    ax.set_xticklabels([])
    ax.set_yticklabels([])

#   Set title to display frame no.
  ax.set_title('Frame no.: '+str(args.frame))
  if args.notitle:
    ax.set_title('')

#   Add a colour bar
  if args.cb:
    plt.colorbar(im)
  
  def on_press(event):
    global ic
    if event.key=='=':
      axes=event.canvas.figure.get_axes()[0]
      if ic != nframes-1:
        ic+=1
      xlim=axes.get_xlim()
      ylim=axes.get_ylim()
      lev_min=shaped_data[ic].min()
      lev_max=shaped_data[ic].max()
      axes.imshow(shaped_data[ic].transpose(),cmap=plt.cm.jet,vmin=lev_min,vmax=lev_max,origin='lower',interpolation='nearest')
      if args.noticks:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
      axes.set_title('Frame no.: '+str(ic+1))
      if args.notitle:
        ax.set_title('')
      axes.set_xlim(xlim)  
      axes.set_ylim(ylim)
      plt.draw()
    if event.key=='-':
      axes=event.canvas.figure.get_axes()[0]
      if ic != 0:
        ic-=1
      xlim=axes.get_xlim()
      ylim=axes.get_ylim()
      lev_min=shaped_data[ic].min()
      lev_max=shaped_data[ic].max()
      axes.imshow(shaped_data[ic].transpose(),cmap=plt.cm.jet,vmin=lev_min,vmax=lev_max,origin='lower',interpolation='nearest')
      if args.noticks:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
      axes.set_title('Frame no.: '+str(ic+1))
      if args.notitle:
        ax.set_title('')
      axes.set_xlim(xlim)  
      axes.set_ylim(ylim)
      plt.draw()


#    print axes[0].get_title()
#    print 'you pressed', event.button, event.xdata, event.ydata
  cid=fig.canvas.mpl_connect('key_press_event',on_press)
#   Close the individual figure
  plt.show()
 

if __name__ == '__main__':
  args = parse_args()
  running(args)



