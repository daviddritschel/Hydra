program kinetic
!  -------------------------------------------------------------------------
!  |   Computes the field of kinetic energy (u^2+v^2)/2 and writes the     |
!  |   output to kk.r4.  Reads data from qq.r4 on the main inversion grid. |
!  -------------------------------------------------------------------------

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use spectral

implicit none

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(ny,nx),uu(ny,nx),vv(ny,nx),pp(ny,nx)

 !Work quantities:
real:: qqr4(ny,nx)
real:: t

integer:: loop,loop1,loop2,iread,igrec,ix,iy

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open file containing PV field:
open(31,file='qq.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

open(41,file='kk.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

!---------------------------------------------------------------
 !Read data and process:

write(*,*) ' Enter range of frames to process (or 0 0 for all):'
read(*,*) loop1,loop2

if (loop1 .eq. 0 .or. loop2 .eq. 0) then
  loop1=1
  loop2=1000000
endif

igrec=0
do loop=loop1,loop2
  iread=0
  read(31,rec=loop,iostat=iread) t,qqr4
  if (iread .ne. 0) exit 
  write(*,'(a,f12.5)') ' Processing t = ',t

   !Define the PV anomaly needed for inversion:
  do ix=1,nx
    do iy=1,ny
      pp(iy,ix)=dble(qqr4(iy,ix))-bety(iy)
    enddo
  enddo

   !Convert to spectral space for inversion:
  call ptospc(nx,ny,pp,qq,xfactors,yfactors,xtrig,ytrig)

   !Invert PV to get velocity field (uu,vv):
  call main_invert(qq,uu,vv,pp)

   !Write diagnostic data for this time:
  do ix=1,nx
    do iy=1,ny
      qqr4(iy,ix)=real(f12*(uu(iy,ix)**2+vv(iy,ix)**2))
    enddo
  enddo
  igrec=igrec+1
  write(41,rec=igrec) t,qqr4
enddo

 !Close files:
close(31)
close(41)

write(*,*)
write(*,*) ' The kinetic energy field is ready in kk.r4'

return
end subroutine

 !End main program
end program
!=======================================================================
