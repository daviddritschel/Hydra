program modes
!  -------------------------------------------------------------------------
!  |   Calculates the streamfunction for each vertical mode from data in   |
!  |   qq1 & qq2.r4.  Writes pm1 & pm2.r4 in the evolution subdirectory.   |
!  -------------------------------------------------------------------------

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use spectral
use generic

implicit none
double precision:: fhb(0:ny,0:nxm1)
double precision:: dum
integer:: ix,iy

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
  open(12,file='topo.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) dum,fhb
  close(12)
else
  do ix=0,nxm1
    do iy=0,ny
      fhb(iy,ix)=zero
    enddo
  enddo
endif

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(0:ny,0:nxm1,nz),pp(0:ny,0:nxm1,nz)
double precision:: uu(0:ny,0:nxm1,nz),vv(0:ny,0:nxm1,nz)
double precision:: p1(0:ny,0:nxm1),p2(0:ny,0:nxm1)

 !Diagnostic quantities:
double precision:: t
real:: q1r4(0:ny,0:nxm1), q2r4(0:ny,0:nxm1), tr4

integer:: loop, iread, ix, iy

!---------------------------------------------------------------
 !Open files containing PV field in each layer:
open(31,file='evolution/qq1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(32,file='evolution/qq2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open files to contain the streamfunction in each mode:
open(41,file='evolution/pm1.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(42,file='evolution/pm2.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,q1r4
  if (iread .ne. 0) exit 
  read(32,rec=loop,iostat=iread) tr4,q2r4

  t=dble(tr4)
  write(*,'(a,f13.5)') ' Processing t = ',t

  do ix=0,nxm1
    do iy=0,ny
      qq(iy,ix,1)=dble(q1r4(iy,ix))
      qq(iy,ix,2)=dble(q2r4(iy,ix))
    enddo
  enddo

   !Invert PV to get velocity field (uu,vv) and streamfunction pp:
  call main_invert(qq,fhb,t,uu,vv,pp)

   !Compute mode 1 and mode 2 parts of psi (in pp):
  do ix=0,nxm1
    do iy=0,ny
      p1(iy,ix)=vec11*pp(iy,ix,1)+vec12*pp(iy,ix,2)
      p2(iy,ix)=vec21*pp(iy,ix,1)+vec22*pp(iy,ix,2)
    enddo
  enddo

   !Write fields at this time:
   !Convert to real*4:
  do ix=0,nxm1
    do iy=0,ny
      q1r4(iy,ix)=real(p1(iy,ix))
      q2r4(iy,ix)=real(p2(iy,ix))
    enddo
  enddo

  write(41,rec=loop) tr4,q1r4
  write(42,rec=loop) tr4,q2r4

enddo

 !Close files:
close(31)
close(32)
close(41)
close(42)

write(*,*)
write(*,*) ' The streamfunctions for mode 1 and 2 are in pm1 & pm2.r4.'

return
end subroutine

 !End main program
end program
!=======================================================================
