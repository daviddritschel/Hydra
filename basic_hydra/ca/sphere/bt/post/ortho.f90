program ortho
!  ---------------------------------------------------------------------
!  |   Creates an orthographic projection of the data in zz.r4 -       |
!  |   writes zzop.r4 which can be imaged using dataview.              |
!  |                                                                   |
!  |   Adapted from r4toc2.f90 on 22 March 2015 by dgd @ St Andrews    |
!  ---------------------------------------------------------------------

use constants

implicit none

real:: rlatc,rlonc,clatc,slatc,b2c2s2,dy
integer:: ny,nz,nbprec
integer:: kbeg,kend,kint

!---------------------------------------------------------------------
write(*,*) ' Latitude & longitude of the direction of view (degrees)?'
read(*,*) rlatc,rlonc

write(*,*) ' Width (in pixels) of output image?'
read(*,*) ny

 !Work out height of output image:
clatc=cos(rlatc*pi/180.)
slatc=sin(rlatc*pi/180.)
b2c2s2=(asp*clatc)**2+slatc**2
dy=2./float(ny)
nz=2*int(0.5+sqrt(b2c2s2)/dy)

 !Define output record length (in bytes):
nbprec=4*(ny*nz+1)

 !Image a time sequence from data in the main job directory:
write(*,*) ' Beginning & ending frames, and the interval (e.g. 1)?'
read(*,*) kbeg,kend,kint

 !Open input data file:
open(11,file='zz.r4',form='unformatted',access='direct', &
        status='old',recl=nbytes)

 !Open output data file:
open(22,file='zzop.r4',form='unformatted',access='direct', &
        status='replace',recl=nbprec)

 !Project data orthographically:
call opgrid(rlatc,rlonc,ny,nz,kbeg,kend,kint)

 !Close files:
close(11)
close(22)

write(*,*)
write(*,*) ' *** Image the results using'
write(*,'(a,i4,1x,i4)') ' dataview zzop.r4 -ndim ',ny,nz

 !Create a file containing the command to image the results:
open(33,file='op_image_command',status='replace')
write(33,'(a,i4,1x,i4)') ' dataview zzop.r4 -ndim ',ny,nz
write(33,*)
write(33,'(a,f5.1,a,f6.1,a)') ' View from ',rlatc, &
                   ' degrees latitude and ',rlonc,' degrees longitude.'
close(33)

!============================================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!============================================================================

subroutine opgrid(rlatc,rlonc,ny,nz,kbeg,kend,kint)
! Projects lat-lon data in zz.r4 orthographically on an ellipsoid

implicit none

integer:: ny,nz,kbeg,kend,kint

integer:: i,j,k,loop
integer:: iy,iz,ic,ip1,jp1

real:: qq(0:ngp1,nt),qqop(nz,ny)
real:: yg(ny),zg(nz)
real:: qqbg,t,dy,dz,b2c2s2,det
real:: rlatc,clatc,slatc,rlonc,clonc,slonc
real:: xp,yp,zp,xm,xt,yt,zt,ri,rj,aa,bb,cc,dd,qqt

!------------------------------------------------------------------
write(*,*) ' Background field level to use for point off the surface?'
read(*,*) qqbg

!------------------------------------------------------------------
 !Define points (yg,zg) in cross-sectional plane of view:
dy=2.0/float(ny)
do iy=1,ny
  yg(iy)=dy*(float(iy)-0.5)-1.0
enddo
dz=dy
do iz=1,nz
  zg(iz)=dz*(float(iz-nz/2)-0.5)
enddo

clatc=cos(rlatc*pi/180.)
slatc=sin(rlatc*pi/180.)

clonc=cos(rlonc*pi/180.)
slonc=sin(rlonc*pi/180.)

b2c2s2=(asp*clatc)**2+slatc**2

!------------------------------------------------------------------
 !Process selected frames:
loop=0
do k=kbeg,kend,kint
   !Read each frame of the data:
  read(11,rec=k) t,qq(1:ng,1:nt)

   !Copy latitudes adjacent to poles (j = 1 and ng) with a pi
   !shift in longitude to simplify interpolation below:
  do i=1,ng
    ic=i+ng
    qq(0,i)=qq(1,ic)
    qq(0,ic)=qq(1,i)
    qq(ngp1,i)=qq(ng,ic)
    qq(ngp1,ic)=qq(ng,i)
  enddo

   !Do interpolation in chosen perspective:
  do iy=1,ny
    yp=yg(iy)
    do iz=1,nz
      zp=zg(iz)
      det=(1.-yp**2)*b2c2s2-zp**2
      if (det .gt. zero) then
        xp=sqrt(det)
        xm=(asp*xp*clatc-zp*slatc)/b2c2s2
         !z/b:
        zt=(asp*zp*clatc+xp*slatc)/b2c2s2
        xm=xp*clatc-zp*slatc
        yt=yp*clonc+xm*slonc
        xt=xm*clonc-yp*slonc

         !Find longitude & latitude then bi-linearly interpolate qq:
        ri=dli*(pi+atan2(yt,xt))
        i=1+int(ri)
        ip1=1+mod(i,nt)
        bb=float(i)-ri
        aa=one-bb

        rj=dli*(hpidl+asin(zt))
        j=int(rj)
        jp1=j+1
        cc=rj-float(j)
        dd=one-cc

        qqop(iz,iy)=bb*(dd*qq(j,i)  +cc*qq(jp1,i))   &
                   +aa*(dd*qq(j,ip1)+cc*qq(jp1,ip1))
      else
        qqop(iz,iy)=qqbg
      endif
    enddo
  enddo

   !Write character image:
  loop=loop+1
  write(22,rec=loop) qqop
enddo

end subroutine

!=========================================================

end program
