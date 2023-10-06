program rossby
! Initialises a Rossby wave

use constants

implicit double precision(a-h,o-z)

double precision,parameter:: k=twopi/ellx,l=twopi/elly
double precision:: qq(ny,nx)

write(*,*) ' The PV is given by q = beta[y+a*sin(nkx+mly)]'
write(*,*) ' where k = 2*pi/L_x and l = 2*pi/L_y.'
write(*,*) ' Enter a:'
read(*,*) a
write(*,*) ' Enter n & m (integer):'
read(*,*) n,m
dn=dble(n)
dm=dble(m)

do ix=1,nx
  x=xmin+glx*dble(ix-1)
  do iy=1,ny
    y=ymin+gly*dble(iy-1)
    qq(iy,ix)=beta*(y+a*sin(dn*k*x+dm*l*y))
  enddo
enddo

open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program
