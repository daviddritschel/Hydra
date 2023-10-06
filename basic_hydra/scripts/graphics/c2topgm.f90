program c2topgm
!===================================================================
 !Converts a selected character image consisting of 2 characters per
 !pixel and renders to a grayscale image in a PGM file format.
!===================================================================

 !Pass arguments in by a script:
integer,parameter:: nx = N_X, ny = N_Y, nbpc=N_BYTES_PER_CHAR
 !npbc: number of bytes used to represent a single character (1 or 4)

integer,parameter:: npix=nx*ny

 !Various character strings:
character(len=2),parameter:: magic='P5'
character(len=3),parameter:: levs ='255'
character(len= 5):: nxstr,nystr
character(len= 2):: v(ny,nx)
character(len= 1):: bmap(nx,ny)
character(len= 9):: otail
character(len=40):: prefix, ifile, ofile

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

write(*,*) ' Frame or record to convert to pgm?'
read(*,*) irec

 !Construct output file name and open it:
ofile = prefix
otail = '_0000.pgm'
write(otail(2:5),'(i4.4)') irec
write(ofile(nc+1:nc+9),'(a9)') otail
write(*,*)

 !Save file name for use by calling script
open(20, file = 'pgmfilename', status = 'unknown')
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

 !Construct grayscale bitmap:
do iy=1,ny
  do ix=1,nx
    ich=int((256*ichar(v(iy,ix)(1:1))+ichar(v(iy,ix)(2:2)))/255)
    bmap(ix,ny+1-iy) = char(ich)
  enddo
enddo

write(11) bmap
close(11)

end program
