program c2tops
!===================================================================
 !Converts a selected character image (having two characters per
 !pixel to render 16 bit colour (256^2 colours)) to postscript.

 !From mad6n@fermi.clas.virginia.edu Thu Jun 15 17:23:13 1995
 !(due originally to Michael A Dritschel; adapted by his brother, 
 ! David, on Wed Nov 16 2005 ... a little bit later.)

 !Final version completed 23 June 2013 by dgd
!===================================================================

 !Pass arguments in by a script:
integer,parameter:: nx = N_X, ny = N_Y, nbpc=N_BYTES_PER_CHAR
 !npbc: number of bytes used to represent a single character (1 or 4)
integer,parameter:: nppi = 72, npix = nx*ny
 !nppi: number of pixels/inch assumed in the postscript translation
real,parameter:: a4width = 8.27, a4height = 11.69
 !Note: A4 is 8.27 x 11.69 inches

 !Colourmap enties:
integer,parameter:: ncol=256*256
integer:: red(0:ncol-1), green(0:ncol-1), blue(0:ncol-1)

 !Various character strings:
character(len= 1),dimension(16),parameter:: hexdig=['0','1','2','3','4', &
                          & '5','6','7','8','9','A','B','C','D','E','F']
character(len= 2):: vc(3*npix),v(npix)
character(len= 8):: otail
character(len=13):: npistr
character(len=27):: scastr
character(len=31):: trastr
character(len=55):: bobstr
character(len=40):: prefix, ifile, ofile
character(len=64):: blanks, vctemp

 !--------------------------------------------------------------
 !Set up basic constants and color/grey maps:
width = float(nx)/float(nppi)

write(*,*) ' Image width (in inches)?'
write(*,'(a,f9.6)') ' ===> Suggested value = ',width
read(*,*) width

yxasp = float(ny)/float(nx)
height = yxasp*width

write(*,*)
write(*,'(a,f10.6,a)') ' The image width will be ',width,' in.'
write(*,'(a,f10.6,a)') '  and its height will be ',height,' in.'
write(*,*)

do i = 1, 64
  blanks(i:i) = ' '
enddo

write(*,*) ' (1) Grey-scale or (2) colour postscript? '
read(*,*) iopt

if (iopt .eq. 1) then
   !Generate linear greymap in the red channel:
  do i = 0, ncol-1
    red(i)=i
  enddo
else
   !Read from existing colourmap:
  open(9,file='colourmap',status='unknown')
  do i = 0, ncol-1
    read(9,*) r,g,b
    red(i)  =int(r*ncol)
    green(i)=int(g*ncol)
    blue(i) =int(b*ncol)
  enddo
  close(9)
endif

 !--------------------------------------------------------------
 !Get input & output file names and open the files:
write(*,*) ' Input file prefix (before .c2)? '
read(*,'(a)') prefix

 !Number of characters in string without trailing spaces:
nc=len(trim(prefix))

ifile = prefix
write(ifile(nc+1:nc+3),'(a3)') '.c2'

open(10,file = ifile,form = 'unformatted',access = 'direct', status = 'old', &
      & recl = 2*npix/nbpc)
 !recl = npix/4 for ifort; = npix for gfortran etc.

write(*,*) ' Frame or record to convert to postscript?'
read(*,*) irps

 !Construct output file name and open it:
ofile = prefix
otail = '_0000.ps'
write(otail(2:5),'(i4.4)') irps
write(ofile(nc+1:nc+8),'(a8)') otail
write(*,*)

open(11,file = ofile,status = 'unknown')

 !--------------------------------------------------------------
 !Read chosen frame/record from input data file:
read(10,rec = irps) v
close(10)

 !--------------------------------------------------------------
 !Convert characters (in v) to hexidecimals (in vc):

 !Here hexdig contains the character representation of the hexidecimal
 !digits, vc is the matrix for the image (represented as one long
 !string: npix is the number of bins in the image, ie, the resolution in
 !the x direction times the resolution in the y direction), ofile is used
 !to store the name of the outputted postscript file, blanks is a line
 !of 64 blanks, and vctemp is temporary storage.

if (iopt .eq. 1) then
   !Grey-scale postscript:
  do i = 1, npix
    ich1 = ichar(v(i)(1:1))
    ich2 = ichar(v(i)(2:2))
    ich=256*ich1+ich2
    icdens = red(ich)
    ihex1  = icdens/4096 + 1
    ihex2 = (icdens-(ihex1-1)*4096)/256 + 1
    vc(i) = hexdig(ihex1)//hexdig(ihex2)
  enddo

   !Work out number of lines to write to ps file:
  nlines = npix/32
  nremain = npix - nlines*32

else
   !Colour postscript:
  do i = 1, npix
    ich1 = ichar(v(i)(1:1))
    ich2 = ichar(v(i)(2:2))
    ich=256*ich1+ich2
     !Red:
    icdens = red(ich)
    ihex1  = icdens/4096 + 1
    ihex2 = (icdens-(ihex1-1)*4096)/256 + 1
    vc(3*i-2) = hexdig(ihex1)//hexdig(ihex2)
     !Green:
    icdens = green(ich)
    ihex1  = icdens/4096 + 1
    ihex2 = (icdens-(ihex1-1)*4096)/256 + 1
    vc(3*i-1) = hexdig(ihex1)//hexdig(ihex2)
     !Blue:
    icdens = blue(ich)
    ihex1  = icdens/4096 + 1
    ihex2 = (icdens-(ihex1-1)*4096)/256 + 1
    vc(3*i) = hexdig(ihex1)//hexdig(ihex2)
  enddo

   !Work out number of lines to write to ps file:
  nlines = 3*npix/32
  nremain = 3*npix - nlines*32

endif

 !In addition one could put coordinate axes and so forth in without much
 !trouble.  So to draw the unit circle, use:

 !    do iang = 1,nangles+1
 !      t = (iang-1)*twopi/float(nangles)
 !      ix = 1 + int(ddfac*(cos(t)+1))
 !      iy = 1 + int(ddfac*(sin(t)+1))
 !      vc(nx*iy + ix) = '00'
 !    enddo

 !--------------------------------------------------------------
 !Now here is how we write the postscript files.
write(*,*) ' Writing the postscript image file, '//ofile
write(*,*)

 !Save file name for use by calling script
open(20, file = 'psfilename', status = 'unknown')
rewind 20
write(20,'(a)') ofile
close(20)

 !Postscipt files have lines with no more than 64 characters across, 
 !so nlines is the number of lines devoted to the data in the postscript 
 !file (with some remaining characters given by nremain if the length of 
 !vc is not divisible by 16 (each entry of vc is four characters).  
 !We write the files sequentially rather than directly since I could only 
 !see how to produce the files using formatted statements.

 !Here is the code for the postscript preamble:

write(11,'(a)') '%!'

bobstr = '%%BoundingBox:                                         '
x1 = (a4width - width)/2.
write(bobstr(17:25),'(f9.4)') float(nppi)*x1
x2 = (a4width + width)/2.
write(bobstr(27:35),'(f9.4)') float(nppi)*x2
y1 = (a4height - height)/2.
write(bobstr(37:45),'(f9.4)') float(nppi)*y1
y2 = (a4height + height)/2.
write(bobstr(47:55),'(f9.4)') float(nppi)*y2
write(11,'(a)') bobstr
write(11,'(a)') '%%EndComments'
write(11,'(a)') 'gsave'

npistr = '        scale'
write(npistr(1:3),'(i3)') nppi
write(npistr(5:7),'(i3)') nppi
write(11,'(a)') npistr

if (iopt .eq. 1) then
  write(11,'(a,i7,a)') '/imline ',2*nx,' string def'
else
  write(11,'(a,i7,a)') '/imline ',6*nx,' string def'
endif

write(11,'(a)') '/drawimage {'
write(11,'(a,i6,a,i6,a)') '    ',nx,' ',ny,' 8'
write(11,'(a,i6,a,i7,a,i6,a)') '     [',nx,' 0 0 ',-ny,' 0 ',ny,']'
write(11,'(a)') '     {currentfile imline readhexstring pop}'

if (iopt .eq. 1) then
  write(11,'(a)') '     image} def'
else
  write(11,'(a)') ' false 3 colorimage} def'
endif

trastr = '                      translate'
write(trastr( 1: 9),'(f9.6)')  x1
write(trastr(11:21),'(f11.6)') a4height - y1
write(11,'(a)') trastr

scastr = '                      scale'
write(scastr( 1: 9),'(f9.6)')   width
write(scastr(11:21),'(f11.6)') -height
write(11,'(a)') scastr

write(11,'(a)') 'drawimage'

 !A few comments on the above.  It is assumed that the data is 8 bit and 
 !that the image is to be read from the top left right and down.  

 !To put the data in the file, we use the following:
do il = 1,nlines
  istep = (il-1)*32
  vctemp = blanks
  do ik = 1,32
    vctemp(2*ik-1:2*ik) = vc(istep+ik)
  enddo
  write(11,'(a64)') vctemp
enddo

if (nremain .ne. 0) then
  istep = (nlines-1)*32
  vctemp = blanks
  do ik = 1,nremain
    vctemp(2*ik-1:2*ik) = vc(istep+ik)
  enddo
  write(11,'(a64)') vctemp
endif

 !This takes successive 64 character pieces of vc, puts them in vctemp,
 !and writes them to the file.  The last conditional takes care of the
 !case when the last segment is not 64 characters long.

 !Finally, the postamble(?) is given by

write(11,'(a)') 'showpage'
write(11,'(a)') 'grestore'

close(11)

end program
