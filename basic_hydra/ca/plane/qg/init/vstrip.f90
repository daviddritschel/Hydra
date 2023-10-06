program strip
! Initialises a PV strip.

use constants

implicit double precision(a-h,o-z)

!integer,parameter:: ngridf=nxf*nyf,nbytesf=4*(ngridf+1)
double precision,parameter:: qmax=4.d0*pi
double precision:: qod0(nyu/2),qod1(nyu/2),qod2(nyu/2)
double precision:: qev0(0:nyu/2),qev1(0:nyu/2),qev2(0:nyu/2)
double precision:: qa(nyu,nxu),qq(ny,nx)

write(*,*) ' Enter the width of the strip:'
read(*,*) wid
hwid=wid/two

write(*,*) ' The upper edge is displaced by A_1*sin(m_1*kx)+A_2*sin(m_2*kx)' 
write(*,*) ' where k = 2*pi/L_x.'
write(*,*) ' Enter A_1 & A_2:'
read(*,*) amp1,amp2
write(*,*) ' Enter m_1 & m_2:'
read(*,*) m1,m2
fm1=dble(m1)
fm2=dble(m2)

glxf=ellx/dble(nxu)
glyf=elly/dble(nyu)
akx=twopi/ellx
do ix=1,nxu
  x=xmin+glxf*dble(ix-1)
  y1=-hwid
  y2=hwid+amp1*sin(fm1*akx*x)+amp2*sin(fm2*akx*x)
  do iy=1,nyu
    y=ymin+glyf*dble(iy-1)
    if ((y2-y)*(y-y1) .gt. zero) then
      qa(iy,ix)=4.d0*qmax*(y2-y)*(y-y1)/(y2-y1)**2
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

open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program
