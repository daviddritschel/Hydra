program energy
!  -------------------------------------------------------------------------
!  |   Computes the integrated kinetic, potential and total energies.      |
!  |                                                                       |
!  |   Output is to the formatted file "kepete.asc" which contains data    |
!  |   of the form                                                         |
!  |           t, <u^2+v^2>/2, <b^2>/2, <u^2+v^2+b^2>/2                    | 
!  |   with each record corresponding to a given time read from the input  |
!  |   file qq.r4.
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
double precision:: qa(ny,nx),qq(ny,nx),uu(ny,nx),vv(ny,nx)
 !Spectral fields:
double precision:: qs(nx,ny),pp(nx,ny)

 !Diagnostic quantities:
real:: qqr4(ny,nx)
real:: t
double precision:: ke,pe

integer:: loop,iread,ix,iy

!---------------------------------------------------------------
 !Open file containing buoyancy field:
open(31,file='qq.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

 !Open output files:
open(11,file='kepete.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) t,qqr4
  if (iread .ne. 0) exit 
  write(*,'(a,f12.5)') ' Processing t = ',t

   !Define the buoyancy anomaly needed for inversion:
  qq=dble(qqr4)
  qa=qq

   !Convert qa to spectral space as qs (note, qa is modified):
  call ptospc(nx,ny,qa,qs,xfactors,yfactors,xtrig,ytrig)

   !Invert buoyancy to get velocity field (uu,vv):
  call main_invert(qs,uu,vv,pp)

   !Compute energy components:
  ke=zero
  pe=zero
  do ix=1,nx
    do iy=1,ny
      ke=ke+uu(iy,ix)**2+vv(iy,ix)**2
      pe=pe+qq(iy,ix)**2
    enddo
  enddo
  ke=f12*garea*ke
  pe=f12*garea*pe

   !Write diagnostic data:
  write(11,'(f12.5,3(1x,f16.11))') t,ke,pe,ke+pe
enddo

 !Close files:
close(11)
close(31)

write(*,*)
write(*,*) ' The results are ready in kepete.asc'

return
end subroutine

 !End main program
end program
!=======================================================================
