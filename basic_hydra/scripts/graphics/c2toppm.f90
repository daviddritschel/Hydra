program c2toppm
!===================================================================
 !Converts a selected character image consisting of 2 characters per
 !pixel and renders using colourmap to a true colour PPM file.
!===================================================================

 !Pass arguments in by a script:
integer,parameter:: nx = N_X, ny = N_Y, nbpc=N_BYTES_PER_CHAR
 !npbc: number of bytes used to represent a single character (1 or 4)

integer,parameter:: npix=nx*ny

 !Colourmap enties:
integer,parameter:: ncol=256*256
integer:: red(0:ncol-1), green(0:ncol-1), blue(0:ncol-1)

 !Various character strings:
character(len=2),parameter:: magic='P6'
character(len=3),parameter:: levs ='255'
character(len= 5):: nxstr,nystr
character(len= 2):: v(ny,nx)
character(len= 3):: bmap(nx,ny)
character(len= 9):: otail
character(len=40):: prefix, ifile, ofile

 !--------------------------------------------------------------
 !Read from existing colourmap:
open(9,file='colourmap',status='unknown')
do i = 0, ncol-1
  read(9,*) r,g,b
  red(i)  =int(r*256.0)
  green(i)=int(g*256.0)
  blue(i) =int(b*256.0)
enddo
close(9)

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

write(*,*) ' Frame or record to convert to ppm?'
read(*,*) irec

 !Construct output file name and open it:
ofile = prefix
otail = '_0000.ppm'
write(otail(2:5),'(i4.4)') irec
write(ofile(nc+1:nc+9),'(a9)') otail
write(*,*)

 !Save file name for use by calling script
open(20, file = 'ppmfilename', status = 'unknown')
write(20,'(a)') ofile
close(20)

open(11,file=ofile,status='unknown',form='unformatted',access='stream')

 !Convert dimensions to character strings:
write(nxstr,'(i5.5)') nx
write(nystr,'(i5.5)') ny

 !Write header:
write(11) magic//' '
write(11) nxstr//' '
write(11) nystr//' '
write(11) levs//' '

 !--------------------------------------------------------------
 !Read chosen frame/record from input data file:
read(10,rec = irec) v
close(10)

 !Construct rgb bitmap:
do iy=1,ny
  do ix=1,nx
    ich=256*ichar(v(iy,ix)(1:1))+ichar(v(iy,ix)(2:2))
    bmap(ix,ny+1-iy) = char(red(ich))//char(green(ich))//char(blue(ich))
  enddo
enddo

write(11) bmap
close(11)

end program
