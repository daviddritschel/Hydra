program vort

use constants

! This routine sets up initial buoyancy and vorticity fields for
! casl in a rectangular aperiodic domain.
! The initial condition is a Gaussian vortex in a stratified domain. 
! Adapted by sek on 5 Sept 2012 from an original f77 code by dgd @ St Andrews
! Adapted further by jm on 19 Jun 2014

implicit double precision (a-h,o-z)

 !Local vorticity and buoyancy arrays:
double precision:: zz(0:ny,0:nx),bb(0:ny,0:nx)

 !Centre of mass
double precision:: x_cen,y_cen

 !Gaussian parameters
double precision:: A,epsilon

write(*,'(a,2(f6.2))') ' We consider a domain of height and length: ',elly,ellx
write(*,'(a,2(i6))') ' The y and x grid resolution is: ',ny,nx
gridrat=dble(nx)*elly/(dble(ny)*ellx)
write(*,'(a,2(f6.2))') ' Note, the y:x grid length ratio is: ',gridrat
write(*,*)

write(*,*) ' Enter the height of the higher density zone:'
read(*,*) h

write(*,*) ' Enter the width of the density transition zone:'
read(*,*) d

write(*,*) ' Enter the x-value of the centre of vorticity:'
read(*,*) x_cen

write(*,*) ' Enter the y-value of the centre of vorticity:'
read(*,*) y_cen

write(*,*) ' Enter the value of the amplitude of the Gaussian vortex:'
read(*,*) A

write(*,*) ' Enter the value of the scale width of the vortex:'
read(*,*) epsilon

 !Set up buoyancy distribution:
do ix=0,nx
  x=glx*dble(ix)
  do iy=0,ny
    y=gly*dble(iy)
    zz(iy,ix)=A*exp(-((x-x_cen)**2+(y-y_cen)**2)/epsilon)
    bb(iy,ix)=f12*erfc((y-h)/d)
  enddo
enddo

 !Write vorticity distribution to file:
open(20,file='zz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

 !Write buoyancy distribution to file:
open(20,file='bb_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,bb
close(20)
      
end program
