program init_tracer
! Initialises the tracer anomaly field as A/(dx*dy) at a specified
! grid point ix,iy.

use constants

implicit none

double precision:: cc(ny,nx)
double precision:: a
integer:: ix,iy

write(*,*)
write(*,*) ' We take the tracer anomaly to have the form of a delta'
write(*,*) ' function spread over one grid box centred on a grid point,'
write(*,*) ' i.e. c(x,y,0) = A*delta(x,y).'
write(*,*)
write(*,*) ' Enter A:'
read(*,*) a

write(*,*) ' Enter the grid point location, ix & iy:'
read(*,*) ix,iy

cc=zero
cc(iy,ix)=a/garea

open(11,file='cc_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,cc
close(11)

end program init_tracer
