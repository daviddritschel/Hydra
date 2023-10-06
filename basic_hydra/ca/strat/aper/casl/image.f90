program mgrid

!  ------------------------------------------------------------------
!  |   Creates images of bouyancy or vorticity from data in bb.dat  |
!  |   or zz.dat.                                                   |
!  ------------------------------------------------------------------

use constants

implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: ix1=0,ix2=nx,iy1=0,iy2=ny
integer,parameter:: npx=ix2+1-ix1,npy=iy2+1-iy1
integer,parameter:: npix=npx*npy

integer,parameter:: nbpc=1
 !  npbc: number of bytes used to represent a single character;
 !        typically 1 when compiled with gfortran but can be 4 

double precision:: qq(0:ny,0:nx)
character bqq(npix)*1
character black*1
character prefix*20
character infile*24
character outfile*24
logical remove

black=char(0)

!----------------------------------------------------------------

write(*,*)
write(*,*) ' Enter the file prefix (before .dat, e.g. qq):'
read(*,'(a)') prefix
nchar=1
do while (prefix(nchar:nchar) .ne. " ") 
  nchar=nchar+1
enddo
nchar=nchar-1

write(infile(1:nchar),'(a)') prefix(1:nchar)
write(infile(nchar+1:nchar+4),'(a)') '.dat'
write(outfile(1:nchar),'(a)') prefix(1:nchar)
write(outfile(nchar+1:nchar+4),'(a)') '.bin'

 !Open input data file:
open(15,file=infile(1:nchar+4),status='old')

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
  qqmax=1.0000001d0*qqmax
  qqmin=1.0000001d0*qqmin
  range=qqmax-qqmin
  qadd=10.d0*range-qqmin
  cfac=240.d0/range
endif

 !Open character movie file(s):
open(21,file=outfile(1:nchar+4),form='unformatted',&
     &access='direct',status='unknown',recl=npix/nbpc)
 !recl = npix/4 for ifort; = npix for gfortran etc.

loop=0
ierr=0
 !Read field and process:
do  
  call readframe(ierr)
  if (ierr .eq. 1) exit 
  q1=1.d14
  q2=-q1
  do ix=0,nx
    do iy=0,ny
      q1=min(q1,qq(iy,ix))
      q2=max(q2,qq(iy,ix))
    enddo
  enddo
  if (jopt .eq. 1) then
    qqmax=1.0000001d0*max(abs(q1),abs(q2))
    qqmin=-qqmax
    range=qqmax-qqmin
    qadd=10.d0*range-qqmin
    cfac=240.d0/range
  else if (jopt .eq. 2) then
    qqdif=q2-q1
    qqmax=q2+0.0000001d0*qqdif
    qqmin=q1-0.0000001d0*qqdif
    range=qqmax-qqmin
    qadd=10.d0*range-qqmin
    cfac=240.d0/range
  endif

  write(*,'(a,f11.5,2(a,f13.7))') ' t = ',t,'    q_min = ',q1,'    q_max = ',q2

 !Write pixel map for this frame:
  loop=loop+1

 !Convert field to character values:
  ixoff=1-ix1
  do ix=ix1,ix2
    do iy=iy1,iy2
      bqq(npx*(iy2-iy)+ixoff+ix)=char(7+nint(cfac*mod(qq(iy,ix)+qadd,range)))
    enddo
  enddo
  write(21,rec=loop) (bqq(ii),ii=1,npix)

enddo
close(21)

write(*,*)
write(*,*) ' Data is ready for display.'
write(*,*)
write(*,*) ' *** Use the command'
write(*,*)
write(*,'(2(a,i4),a)') ' xvidi -dwidth ',npx,' -dheight ',npy,&
                      &' '//outfile(1:nchar+4)//' &'
write(*,*)


contains 

!=========================================================

subroutine readframe(ierr)

! Subroutine reads a frame of the file. Exits with ierr=0 if 
! the frame was read successfully, 1 otherwise.

implicit double precision(a-h,o-z)
implicit integer(i-n)

ierr=0
read(15,*,iostat=itmp) t
if (itmp .ne. 0) then
  ierr=1
  return
endif

do ix=0,nx
  do iy=0,ny
    read(15,*,iostat=itmp) qq(iy,ix)
    if (itmp .ne. 0) then 
      ierr=1
      return 
    endif
  enddo
enddo

end subroutine

!=========================================================

end program
