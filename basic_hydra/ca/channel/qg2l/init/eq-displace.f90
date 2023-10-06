program displace
! --------------------------------------------------------------------
! |   This routine sets up the thermal equilibrium interface         |
! |   displacements corresponding to the initial PV distribution.    |
! |   This is used to relax the PV back to its initial distribution. |
! |   Only used if the logical "heating" is true.                    |
!---------------------------------------------------------------------

 !Import constants, parameters and spectral module:
use constants
use spectral

implicit none

double precision:: dd(0:ny,0:nxm1)
double precision:: qq(0:ny,0:nxm1,nz),pp(0:ny,0:nxm1,nz)
double precision:: uu(0:ny,0:nxm1,nz),vv(0:ny,0:nxm1,nz)
double precision:: t
integer:: iz

 !----------------------------------------------------------
 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
  open(12,file='topo.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) t,dd
  close(12)
else
  dd=zero
endif

 !----------------------------------------------------------
 !Read in initial PV:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
do iz=1,nz
  read(11,rec=iz) t,qq(:,:,iz)
enddo
close(11)

 !----------------------------------------------------------
 !Initialise spectral module and compute initial flow:
call init_spectral

 !Obtain streamfunction pp in both layers:
call main_invert(qq,dd,t,uu,vv,pp)

 !Store layer 1 interface displacement in dd for writing below:
dd=h1h2kdbarsq*(pp(:,:,1)-alpha*pp(:,:,2))

open(12,file='disp1eq.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(12,rec=1) zero,dd
close(12)

 !If no barotropic mode is present, also write layer 2 displacement:
if (.not. barot) then
  dd=h1h2ackdbarsq*pp(:,:,2)

  open(12,file='disp2eq.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(12,rec=1) zero,dd
  close(12)
endif

end program
