program ellipse_topo
! ---------------------------------------------------------------------|
! |   This routine sets up an elliptically-shaped bottom topography    |
!----------------------------------------------------------------------|

 !Import constants & parameters:
use constants

implicit none

double precision:: bb(0:ny,0:nxm1)
double precision:: amp,ael,bel,phi,cop,sip
double precision:: xg,yg,xx,yy,bavg
integer:: ix,iy

 !Add topography (NOTE: this must have zero domain mean):
write(*,*) ' We consider a topographic PV of the form'
write(*,*) '   q_b = A*exp(-(X/a)^2-(Y/b)^2)'
write(*,*) ' where X & Y are rotated CCW by an angle phi from the x & y axes,'
write(*,*) ' taken to cross through the domain centre.'
write(*,*)
write(*,*) ' Enter A/pi:'
read(*,*) amp
amp=pi*amp
write(*,*) ' Enter a:'
read(*,*) ael
write(*,*) ' Enter b:'
read(*,*) bel
write(*,*) ' Enter phi:'
read(*,*) phi
phi=phi*pi/180.d0

 !---------------------------------------------------------------------------
 !Generate and write topographic PV (NOTE: this must have zero domain mean):
 !NOTE: the centre of the channel is y = 0:
cop=cos(phi)
sip=sin(phi)
do ix=0,nxm1
  xg=xmin+glx*dble(ix)
  do iy=0,ny
    yg=ymin+gly*dble(iy)
    xx=xg*cop+yg*sip
    yy=yg*cop-xg*sip
    bb(iy,ix)=amp*exp(-(xx/ael)**2-(yy/bel)**2)
  enddo
enddo

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

open(11,file='topo.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,bb
close(11)

end program ellipse_topo
