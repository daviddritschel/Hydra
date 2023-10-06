program staircase
! Initialises a PV staircase + a vortex

use constants

implicit double precision(a-h,o-z)

!integer,parameter:: ngridf=nxf*nyf,nbytesf=4*(ngridf+1)
double precision,parameter:: qmax=4.d0*pi
double precision:: qod0(nyu/2),qod1(nyu/2),qod2(nyu/2)
double precision:: qev0(0:nyu/2),qev1(0:nyu/2),qev2(0:nyu/2)
double precision:: qa(nyu,nxu),qq(ny,nx)

write(*,*) ' The background PV is a staircase with a jet spacing of'
write(*,*) ' beta*L_y/2 (2 steps).  A vortex is placed at x = 0 and'
write(*,*) ' y = y_0, having radius R_0 and uniform PV q_0.'
write(*,*) ' Enter y_0, R_0 and 2*q_0/beta*L_y:'
read(*,*) y0,r0,q0nd
q0=beta*hly*q0nd
r0sq=r0**2

b=elly/4.d0
qbg=beta*elly/two

 !Generate PV anomaly, q - beta*y:
glxf=ellx/dble(nxu)
glyf=elly/dble(nyu)
do ix=1,nxu
  x=xmin+glxf*dble(ix-1)
  xsq=x**2
  do iy=1,nyu
    y=ymin+glyf*dble(iy-1)
     !background PV staircase:
    if (y .lt. -b) then
      qa(iy,ix)=-qbg-beta*y
    else if (y .gt. b) then
      qa(iy,ix)= qbg-beta*y
    else
      qa(iy,ix)=-beta*y
    endif
     !Vortex:
    if (xsq+(y-y0)**2 .le. r0sq) then
      qa(iy,ix)=qa(iy,ix)+q0
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

 !Add back beta*y:
do ix=1,nx
  do iy=1,ny
    qq(iy,ix)=qa(iy,ix)-qavg+beta*(ymin+gly*dble(iy-1))
  enddo
enddo

open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program
