program bcjet

!========================================================================
! This routine sets up the thermal equilibrium interface displacement 
! (multiplied by f_0/(H1+H2)) between the two layers, for a baroclinic
! jet (as in jet.f90) together with a state of rest.
! Written 19 Oct 2013 by dgd @ St Andrews
!
!       *** WARNING: Assumes there is a barotropic mode ***
!========================================================================

use constants
use spectral

implicit double precision (a-h,o-z)

 !Local arrays:
double precision:: fhb(0:ny,0:nxm1),qbc(0:ny,0:nxm1)
double precision:: qq(0:ny,0:nxm1,nz),pp(0:ny,0:nxm1,nz)
double precision:: uu(0:ny,0:nxm1,nz),vv(0:ny,0:nxm1,nz)

if (.not. barot) then
  write(*,*) ' *** init/thermal-bcjet.f90 assumes there is a barotropic mode ***'
  write(*,*) ' *** rerun choosing c_rho = 0 ***'
  write(*,*) ' *** stopping! ***'
  stop
endif

!---------------------------------------------------------
 !Initialise inversion constants and arrays:

 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
  open(12,file='topo.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) dum,fhb
  close(12)
else
  do ix=0,nxm1
    do iy=0,ny
      fhb(iy,ix)=zero
    enddo
  enddo
endif

call init_spectral

write(*,'(a,2(f6.2))') ' We consider a domain of height and length: ',elly,ellx
write(*,'(a,2(i6))') ' The y and x grid resolution is: ',ny,nx
gridrat=dble(nx)*elly/(dble(ny)*ellx)
write(*,'(a,2(f6.2))') ' Note, the y:x grid length ratio is: ',gridrat
write(*,*)

write(*,*) ' Eastward or westward jet (1 or -1)?'
read(*,*) idir
 !Choose dq = +/-4*pi:
dq=sign(one,dble(idir))*four*pi
hdq=dq/two

write(*,*) ' Width of the jet?'
read(*,*) w
hw=w/two

 !BC PV gradient in jet:
gamma=dq/w

write(*,*) ' We perturb the centreline of the jet by eps*sin(2*pi*x/L_x).'
write(*,*) ' Enter eps:'
read(*,*) eps

 !Set up the baroclinic PV distribution:
do ix=0,nxm1
  x=xmin+glx*dble(ix)
  dy=eps*sin(twopi*x/ellx)
  y1=-hw+dy
  y2= hw+dy
  do iy=0,ny
    y=ymin+gly*dble(iy)
    if (y .lt. y1) then
      qbc(iy,ix)=-hdq
    else if (y .gt. y2) then
      qbc(iy,ix)= hdq
    else
      qbc(iy,ix)=gamma*(y-dy)
    endif
  enddo
enddo

 !Convert to PV in each layer (add on beta*y):
do ix=0,nxm1
  do iy=0,ny
    qq(iy,ix,1)=bety(iy)-h2*qbc(iy,ix)
    qq(iy,ix,2)=bety(iy)+h1*qbc(iy,ix)
  enddo
enddo

 !Invert PV to get velocity field (uu,vv) and streamfunction pp:
call main_invert(qq,fhb,uu,vv,pp)

 !Store interface displacement in qbc for writing below:
do ix=0,nxm1
  do iy=0,ny
    qbc(iy,ix)=h1h2kdbarsq*(pp(iy,ix,1)-pp(iy,ix,2))
  enddo
enddo

 !Write equilibrium interface displacement:
open(12,file='disp1eq.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(12,rec=1) dum,qbc
close(12)

 !Write initial undisturbed PV distribution:
do ix=0,nxm1
  do iy=0,ny
    qq(iy,ix,1)=bety(iy)
    qq(iy,ix,2)=bety(iy)
  enddo
enddo
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
do iz=1,nz
  write(11,rec=iz) t,qq(:,:,iz)
enddo
close(11)

 !Write input data to file:
open(12,file='input_for_bcjet',status='unknown')
if (idir .gt. 0) then
  write(12,*) '  Eastward jet with dq = 4*pi'
  write(12,'(4x,f6.4,a)')  h2*gamma/beta,' ! h_2*dq/beta*w'
else
  write(12,*) '  Westward jet with dq = -4*pi'
  write(12,'(4x,f6.4,a)') -h1*gamma/beta,' ! -h_1*dq/beta*w'
endif
write(12,'(4x,f9.5,a)')   h1,' ! lower layer fractional thickness'
write(12,'(4x,f9.5,a)')   h2,' ! upper layer fractional thickness'
write(12,'(4x,f9.5,a)') beta,' ! planetary vorticity gradient, beta'
write(12,'(4x,f9.5,a)')    w,' ! width of the jet'
write(12,'(4x,f9.5,a)')  eps,' ! amplitude of centreline displacement'
close(12)
      
end program
