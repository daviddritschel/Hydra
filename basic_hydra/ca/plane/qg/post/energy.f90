program energy
!  -------------------------------------------------------------------------
!  |   Computes the integrated kinetic, potential and total energies.      |
!  |                                                                       |
!  |   Output is to the formatted file "kepete.asc" which contains data    |
!  |   of the form                                                         |
!  |        t, <u^2+v^2>/2, k_d^2<psi^2>/2, <u^2+v^2+k_d^2*psi^2>/2        | 
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
double precision:: qa(ny,nx),qq(ny,nx),uu(ny,nx),vv(ny,nx),zz(ny,nx)
 !Spectral fields:
double precision:: qs(nx,ny),pp(nx,ny),ss(nx,ny)

 !Diagnostic quantities:
real:: qqr4(ny,nx)
real:: t
double precision:: ke,pe

integer:: loop,iread,ix,iy

!---------------------------------------------------------------
 !Open file containing PV field:
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

   !Define the PV anomaly needed for inversion:
  qq=dble(qqr4)
  qa=qq

   !Convert qa to spectral space as qs (note, qa is modified):
  call ptospc(nx,ny,qa,qs,xfactors,yfactors,xtrig,ytrig)

   !Invert PV to get velocity field (uu,vv):
  call main_invert(qs,uu,vv,pp)

   !Compute kinetic energy:
  zz=uu**2+vv**2
  ke=f12*garea*sum(zz)

   !Compute potential energy:
  if (btropic) then
     !kd = 0 and hence there is no potential energy:
    pe=zero
  else
     !Get streamfunction pp in physical space as zz:
    ss=pp
    call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)
    zz=zz**2
    pe=f12*garea*kdsq*sum(zz)
  endif

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
