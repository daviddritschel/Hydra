program random
!----------------------------------------------------------------------
!    Generates a random phased PV anomaly and/or magnetic potential
!    distribution with a spectrum c * k^{2p-3} * exp[-2*(k/k_0)^2].
!    Here, k_0 is the spectral centroid and p (integer) > 1.
!----------------------------------------------------------------------

use spectral

implicit none

double precision:: qq(ny,nx),aa(ny,nx),ff(ny,nx)
double precision:: uu(ny,nx),vv(ny,nx)
double precision:: ss(nx,ny),pp(nx,ny)
double precision:: urms,brms,pq,pa,fac
integer, dimension(:), allocatable :: seed
integer:: kq,ka,k,i,ngen,ix

!-------------------------------------------------------------
! Initialise inversion constants and arrays:
call init_spectral

!-------------------------------------------------------------
write(*,*) ' We assume a spectrum of the form c * k^{2p-1} * exp[-2*(k/k_0)^2]'
write(*,*) ' for both the PV anomaly q and the magnetic potential A.'

write(*,*)
write(*,*) ' For q, enter the rms velocity, u_rms:'
read(*,*) urms
write(*,*) ' Enter p:'
read(*,*) pq
write(*,*) ' Enter k_0:'
read(*,*) kq

write(*,*)
write(*,*) ' For A, enter the rms magnetic field, b_rms:'
read(*,*) brms
write(*,*) ' Enter p:'
read(*,*) pa
write(*,*) ' Enter k_0:'
read(*,*) ka

write(*,*)
write(*,*) ' Enter an integer seed for the random number generator:'
read(*,*) ngen
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=ngen
do i=1,ngen
  call random_seed(put=seed)
enddo

! Generate q:
if (urms .gt. zero) then
  call ranspec(qq,one,pq,kq)
else
  qq=zero
endif

! Normalise to have correct rms velocity:
ff=qq
call ptospc(nx,ny,ff,ss,xfactors,yfactors,xtrig,ytrig)
call main_invert(ss,uu,vv,pp)
fac=urms/sqrt(dsumi*sum(uu**2+vv**2))   !dsumi = 1/dble(nx*ny)
qq=fac*qq

if (beffect) then
   !Add beta*y (in bety below) to define the full PV:
  do ix=1,nx
    qq(:,ix)=qq(:,ix)+bety
  enddo
endif

! Generate A:
if (brms .gt. zero) then
  call ranspec(aa,one,pa,ka)
else
  aa=zero
endif

! Normalise to have correct rms magnetic field:
ff=aa
call ptospc(nx,ny,ff,ss,xfactors,yfactors,xtrig,ytrig)
call gradient(ss,uu,vv)
fac=brms/sqrt(dsumi*sum(uu**2+vv**2))   !dsumi = 1/dble(nx*ny)
aa=fac*aa

! Write data:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

open(11,file='aa_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,aa
close(11)

end program random
