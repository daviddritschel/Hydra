program r4toc2

!  ------------------------------------------------------------------
!  |   Creates double character files of bouyancy or vorticity      |
!  |   data from bb.r4 or zz.r4                                     |
!  ------------------------------------------------------------------

use constants

implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: ix1=0,ix2=nx,iy1=0,iy2=ny
integer,parameter:: npx=ix2+1-ix1,npy=iy2+1-iy1
integer,parameter:: npix=npx*npy,nbread=(nx+1)*(ny+1)+1

real,parameter:: cfrac=1024./(256.*256.)
 !This parameter is used to access the fraction of the colourmap 
 !extending from cfrac to 1-cfrac (to avoid oversaturation)

integer,parameter:: nbpc=1
 !  npbc: number of bytes used to represent a single character;
 !        typically 1 when compiled with gfortran but can be 4 

real:: crange,coff,qqmax,qqmin,range,q1,q2,qqdif
real:: t,qq(0:ny,0:nx)
character(len=2):: bqq(0:ny,0:nx)
character(len=24):: prefix,infile,outfile

!----------------------------------------------------------------
 !Quantities required for converting data into characters:
crange=float(256**2)
coff=cfrac*crange

write(*,*)
write(*,*) ' Enter the file prefix (before .dat, e.g. qq):'
read(*,'(a)') prefix
nchar=1
do while (prefix(nchar:nchar) .ne. " ") 
  nchar=nchar+1
enddo
nchar=nchar-1

write(infile(1:nchar),'(a)') prefix(1:nchar)
write(infile(nchar+1:nchar+3),'(a)') '.r4'
write(outfile(1:nchar),'(a)') prefix(1:nchar)
write(outfile(nchar+1:nchar+3),'(a)') '.c2'

 !Open input data file:
open(15,file=infile(1:nchar+4),form='unformatted', & 
        & access='direct',status='old',recl=nbread*4)

write(*,*)
write(*,*) ' Let "f" stand for the field shown.'
write(*,*)

write(*,*) ' Choose on of the following options:'
write(*,*)
write(*,*) ' (1) image between -|f|_max and |f|_max,'
write(*,*) ' (2) image between f_min and f_max, or'
write(*,*) ' (3) specify min and max values?'
read(*,*) jopt

if (jopt .eq. 3) then
  write(*,'(a,i1,a)') '  Min & max f to image? '
  read(*,*) qqmin,qqmax
  qqmax=1.0000001e0*qqmax
  qqmin=1.0000001e0*qqmin
  range=qqmax-qqmin
  cfac=(crange-2.*coff)/range
endif

 !Open character movie file(s):
open(21,file=outfile(1:nchar+4),form='unformatted',&
     &access='direct',status='replace',recl=2*npix*nbpc)
 !recl = npix*4 for ifort; = npix for gfortran etc.

loop=0
iread=0
 !Read field and process:
do 
  loop=loop+1
  iread=0
  read(15,rec=loop,iostat=iread) t,qq
  if (iread .ne. 0) exit 
  q1=1.d14
  q2=-q1
  do ix=0,nx
    do iy=0,ny
      q1=min(q1,qq(iy,ix))
      q2=max(q2,qq(iy,ix))
    enddo
  enddo
  if (jopt .eq. 1) then
    qqmax=1.0000001e0*max(abs(q1),abs(q2))
    qqmin=-qqmax
    range=qqmax-qqmin
    cfac=(crange-2.*coff)/range
  else if (jopt .eq. 2) then
    qqdif=q2-q1
    qqmax=q2+0.0000001e0*qqdif
    qqmin=q1-0.0000001e0*qqdif
    range=qqmax-qqmin
    cfac=(crange-2.*coff)/range
  endif

  write(*,'(a,f11.5,2(a,f13.7))') ' t = ',t,'    q_min = ',q1,'    q_max = ',q2

   !Convert each data value to a pair of characters:
  do ix=ix1,ix2
    do iy=iy1,iy2
      ich=nint(coff+cfac*(min(qqmax,max(qqmin,qq(iy,ix)))-qqmin))
      ich1=ich/256+1
      ich2=ich-(ich1-1)*256+1
      bqq(iy,ix)=char(ich1)//char(ich2)
    enddo
  enddo
   !Write character image:
  write(21,rec=loop) bqq(iy1:iy2,ix1:ix2)

enddo
close(15)
close(21)

write(*,*)
write(*,*) ' Data is ready for display.'
write(*,*)
write(*,'(2(a,i5))') ' Note, the image has dimensions ',npx,' by ',npy
write(*,*)

end program
