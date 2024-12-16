program eddy
! Initialises an elliptical eddy.

use constants

implicit double precision(a-h,o-z)

!integer,parameter:: ngridf=nxf*nyf,nbytesf=4*(ngridf+1)
double precision:: qod0(nyu/2),qod1(nyu/2),qod2(nyu/2)
double precision:: qev0(0:nyu/2),qev1(0:nyu/2),qev2(0:nyu/2)
double precision:: qa(nyu,nxu),qq(ny,nx)

write(*,*)
write(*,*) ' We take b_0/N = (1 - s)^p where s = (x/x_0)^2 + (y/y_0)^2'
write(*,*) ' for s < 1, and b_0/N = 0 otherwise. *** Use p = 0 to instead'
write(*,*) ' take b_0/N = e^{-s}.'
write(*,*)
write(*,*) ' Enter p, x_0 and y_0:'
read(*,*) pow,x0,y0

xfac=one/x0
yfac=one/y0

glxf=ellx/dble(nxu)
glyf=elly/dble(nyu)
if (pow > zero) then
   do ix=1,nxu
      x=xfac*(xmin+glxf*dble(ix-1))
      do iy=1,nyu
         y=yfac*(ymin+glyf*dble(iy-1))
         ssq=x**2+y**2
         if (ssq < one) then
            qa(iy,ix)=(one-ssq)**pow
         else
            qa(iy,ix)=zero
         endif
      enddo
   enddo
else
   do ix=1,nxu
      x=xfac*(xmin+glxf*dble(ix-1))
      do iy=1,nyu
         y=yfac*(ymin+glyf*dble(iy-1))
         qa(iy,ix)=exp(-x**2-y**2)
      enddo
   enddo
endif

!------------------------------------------------------------------------
 !Average the buoyancy field in qa to the coarser grid (ny,nx):
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
                   +0.125d0*(qev1(iy)+qev1(iy-1)+qod0(iy)+qod2(iy)) &
                   +0.25d0*qod1(iy)
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

 !Save average for plotting purposes:
open(44,file='average_qq.asc',status='replace')
write(44,*) qavg
close(44)

do ix=1,nx
   do iy=1,ny
      qq(iy,ix)=qa(iy,ix)-qavg
   enddo
enddo

open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program eddy
