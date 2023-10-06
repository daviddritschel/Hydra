program narrow
!-----------------------------------------------------------------
!    Generates a narrow-band spectrum PV distribution having a 
!    maximum absolute PV anomaly of 4*pi.
!-----------------------------------------------------------------

use spectral

implicit double precision(a-h,o-z)
!implicit none

! Set the max abs QG PV anomaly = 4*pi:'
double precision, parameter:: qeddy=4.d0*pi

! Local variables:
double precision:: qa(ny,nx)
!double precision:: ens,qamin,qamax,fmult
!integer:: ix,iy

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Initialise narrow-band spectral forcing to use as initial conditions:
call init_forcing

 !Obtain random spectral field in the array dqdt:
call nb_forcing(-1.d0)

! Transform to physical space:
call spctop(nx,ny,dqdt,qa,xfactors,yfactors,xtrig,ytrig)

! Work out max/min values and total pot. enstrophy:
ens=zero
qamin=qa(1,1)
qamax=qamin
do ix=1,nx
  do iy=1,ny
    ens=ens+qa(iy,ix)**2
    qamin=min(qamin,qa(iy,ix))
    qamax=max(qamax,qa(iy,ix))
  enddo
enddo
ens=ens/(two*dble(nx*ny))

! Renormalise PV:
fmult=qeddy/max(abs(qamax),abs(qamin))
do ix=1,nx
  do iy=1,ny
    qa(iy,ix)=fmult*qa(iy,ix)
  enddo
enddo

! Work out max/min values and total pot. enstrophy:
ens=ens*fmult
qamin=qamin*fmult
qamax=qamax*fmult

! Write data:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qa
close(11)

write(*,'(a,f12.5)') ' rms PV = ',sqrt(2.d0*ens)
write(*,'(a,f12.7,a,f11.7)') ' min PV = ',qamin,'  &  max PV = ',qamax

end program narrow
