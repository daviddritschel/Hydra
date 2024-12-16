program init_tracer
! Initialises the tracer anomaly field as A*sin(2*pi*x/L_x)*sin(2*pi*y/L_y)

use constants

implicit none

double precision:: cc(ny,nx)
double precision:: a,x,y,kx,ky
integer:: ix,iy

write(*,*)
write(*,*) ' We take the tracer anomaly to have the form'
write(*,*) ' c = A*sin(2*pi*x/L_x)*sin(2*pi*y/L_y).'
write(*,*)
write(*,*) ' Enter A:'
read(*,*) a

kx=twopi/ellx
ky=twopi/elly

do ix=1,nx
  x=xmin+glx*dble(ix-1)
  do iy=1,ny
    y=ymin+gly*dble(iy-1)
    cc(iy,ix)=a*sin(kx*x)*sin(ky*y)
  enddo
enddo

open(11,file='cc_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,cc
close(11)

end program init_tracer
