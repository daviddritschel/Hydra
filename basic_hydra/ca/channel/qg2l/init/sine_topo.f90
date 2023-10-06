program cosine_topo
! |-------------------------------------------------------|
! |   This routine sets up sinusoidal bottom topography   |
! |-------------------------------------------------------|

 !Import constants & parameters:
use constants

implicit none

double precision:: bb(0:ny,0:nxm1)
double precision:: topamp,k_topo,l_topo,y0,displa,bavg
integer:: ix,iy

 !Add topography (NOTE: this must have zero domain mean):
write(*,*) ' We consider topography H_b of the form '
write(*,*) ' f_0*H_b/H_1 = A*cos(2*pi*l*(y-y_c(x))/L_y),'
write(*,*) ' where H_1 is the depth of the lowest layer and'
write(*,*) ' y_c = y_0 + D*sin(2*pi*k*x/L_x).'
write(*,*)

write(*,*) ' Enter A:'
read(*,*) topamp

write(*,*) ' Enter k & l (integer):'
read(*,*) k_topo,l_topo
k_topo=dble(nint(k_topo))*twopi/ellx
l_topo=dble(nint(l_topo))*twopi/elly

write(*,*) ' Enter y_0/L_y and D/L_y:'
read(*,*) y0,displa
y0=elly*y0
displa=elly*displa

 !Set up topography:
do ix=0,nxm1
  yc=y0+displa*sin(k_topo*(xmin+glx*dble(ix)))
  do iy=0,ny
    bb(iy,ix)=topamp*cos(l_topo*(ymin+gly*dble(iy)-yc))
  enddo
enddo

 !Compute and remove domain mean:
bavg=zero
do ix=0,nxm1
  bavg=bavg+bb(0,ix)+bb(ny,ix)
enddo
bavg=f12*bavg

do ix=0,nxm1
  do iy=1,nym1
    bavg=bavg+bb(iy,ix)
  enddo
enddo

bavg=bavg/dble(nx*ny)
bb=bb-bavg

 !Write topography to a file:
open(11,file='topo.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,bb
close(11)

end program cosine_topo
