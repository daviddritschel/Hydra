program strip
! Initialises a PV strip.

use constants

implicit double precision(a-h,o-z)

double precision,parameter:: qmax=4.d0*pi
double precision:: qod0(nyu/2),qod1(nyu/2),qod2(nyu/2)
double precision:: qev0(0:nyu/2),qev1(0:nyu/2),qev2(0:nyu/2)
double precision:: qa(nyu,nxu),qq(ny,nx)

write(*,*) ' Enter the width of the strip:'
read(*,*) wid
hwid=wid/two

write(*,*) ' The upper edge is displaced by A_2*sin(2kx)+A_3*sin(3kx) where'
write(*,*) ' k = 2*pi/L_x.'
write(*,*) ' Enter A_2 & A_3:'
read(*,*) amp2,amp3

glxf=ellx/dble(nxu)
glyf=elly/dble(nyu)
akx=twopi/ellx
do ix=1,nxu
  x=xmin+glxf*dble(ix-1)
  y1=-hwid
  y2=hwid+amp2*sin(two*akx*x)+amp3*sin(three*akx*x)
  do iy=1,nyu
    y=ymin+glyf*dble(iy-1)
    if ((y2-y)*(y-y1) .gt. zero) then
      qa(iy,ix)=qmax
    else 
      qa(iy,ix)=zero
    endif
  enddo
enddo

!------------------------------------------------------------------------
 !Average the PV field in qa to the coarser grid (ny,nx):
nxh=nxu
nyh=nyu
do while (nxh .gt. nx)
  nxuf=nxh
  nxh=nxh/2
  nyuf=nyh
  nyh=nyh/2
   !Perform nine-point averaging:
  do iy=1,nyh
    miy=2*iy
    qod2(iy)=qa(miy-1,nyuf)
    qev2(iy)=qa(miy,nyuf)
  enddo
  qev2(0)=qa(nyuf,nyuf)
  do ix=1,nxh
    mix=2*ix
    mixm1=mix-1
    do iy=1,nyh
      miy=2*iy
      qod1(iy)=qa(miy-1,mixm1)
      qod0(iy)=qa(miy-1,mix)
      qev1(iy)=qa(miy,mixm1)
      qev0(iy)=qa(miy,mix)
    enddo
    qev1(0)=qev1(nyh)
    qev0(0)=qev0(nyh)
    do iy=1,nyh
      qa(iy,ix)=0.0625d0*(qev0(iy)+qev0(iy-1)+qev2(iy)+qev2(iy-1)) &
              & +0.125d0*(qev1(iy)+qev1(iy-1)+qod0(iy)+qod2(iy)) &
              &   +0.25d0*qod1(iy)
    enddo
    do iy=1,nyh
      qod2(iy)=qod0(iy)
      qev2(iy)=qev0(iy)
    enddo
    qev2(0)=qev0(0)
  enddo
enddo

 !Calculate and remove average qq:
qavg=zero
do ix=1,nx
  do iy=1,ny
    qavg=qavg+qa(iy,ix)
  enddo
enddo
qavg=qavg/dble(nx*ny)

do ix=1,nx
  do iy=1,ny
    qq(iy,ix)=qa(iy,ix)-qavg
  enddo
enddo

 !Write data:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

 !Write zero magnetic potential:
do ix=1,nx
  do iy=1,ny
    qq(iy,ix)=zero
  enddo
enddo
open(12,file='aa_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(12,rec=1) zero,qq
close(12)

end program
