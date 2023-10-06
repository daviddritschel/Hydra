program ortho
!  ---------------------------------------------------------------------
!  |   Creates an orthographic projection of the data in ff.r4 where   |
!  |   ff is any selected field in the evolution subdirectory.         |
!  |   Writes ffop.r4 which can be imaged using odv.                   |
!  |                                                                   |
!  |   Adapted from bt code on => 12.02.20.21 <= by dgd @ St Andrews   |
!  ---------------------------------------------------------------------

use constants

implicit none

integer,parameter:: nfield=6
real:: rlatc,rlonc
integer:: ny,nz,nbprec
integer:: iopt,jopt
character(len=2):: dname(nfield),field(nfield)
character(len=6):: infile
character(len=8):: outfile
character(len=5):: viewfile

data dname /'qq','qq','zz','dd','gg','hh'/
data field /'qq','qa','zz','dd','gg','hh'/

!---------------------------------------------------------------------
write(*,*) ' Which field do you wish to project?'
write(*,*) '   (1)  PV (q),'
write(*,*) '   (2)  PV anomaly (q-f),'
write(*,*) '   (3)  relative vorticity (zeta),'
write(*,*) '   (4)  velocity divergence (delta),'
write(*,*) '   (5)  acceleration divergence (gamma), or'
write(*,*) '   (6)  dimensionless height anomaly (h)?'
write(*,*) ' Option?'
read(*,*) iopt
if (iopt*(nfield+1-iopt) .le. 0) then
  write(*,*) ' *** Not a valid option; exiting!'
  stop
endif
write(*,*) ' (1) Full or (2) balanced fields?'
read(*,*) jopt
if (jopt*(3-jopt) .le. 0) then
  write(*,*) ' *** Not a valid option; exiting!'
  stop
endif

if (jopt .eq. 1) then
  infile=dname(iopt)//'.r4'
  outfile=field(iopt)//'op.r4'
  viewfile=field(iopt)//'op'
else
  infile='b'//dname(iopt)//'.r4'
  outfile='b'//field(iopt)//'op.r4'
  viewfile='b'//field(iopt)//'op'
endif

write(*,*) ' Latitude & longitude of the direction of view (degrees)?'
read(*,*) rlatc,rlonc

!write(*,*) ' Width (in pixels) of output image?'
!read(*,*) ny
ny=ng !need this to be compatible with the odv viewing script [could use ny=nt]
nz=ny

 !Define output record length (in bytes):
nbprec=4*(ny*nz+1)

 !Open input data file:
open(11,file='evolution/'//infile, form='unformatted',access='direct', &
                                 status='old',recl=nbytes)

 !Open output data file:
open(22,file='evolution/'//outfile,form='unformatted',access='direct', &
                                 status='replace',recl=nbprec)

 !Project data orthographically:
call opgrid(rlatc,rlonc,ny,nz,iopt)

 !Close files:
close(11)
close(22)

write(*,*)
write(*,*) ' *** View the results using'
write(*,*)
write(*,*) ' odv '//viewfile
write(*,*)

!============================================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!============================================================================

subroutine opgrid(rlatc,rlonc,ny,nz,iopt)
! Projects lat-lon data in zz.r4 orthographically on an ellipsoid

implicit none

integer:: ny,nz

integer:: i,j,loop,iread,iopt
integer:: iy,iz,ic,ip1,jp1

real:: qq(0:ngp1,nt),qqop(nz,ny)
real:: yg(ny),zg(nz),cof(ng)
real:: qqbg,t,dy,dz,det
real:: rlatc,clatc,slatc,rlonc,clonc,slonc
real:: xp,yp,zp,xm,xt,yt,zt,ri,rj,aa,bb,cc,dd,qqt

!------------------------------------------------------------------
write(*,*) ' Background field level (for points off the surface)?'
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

if (iopt .eq. 2) then
  do j=1,ng
    cof(j)=fpole*sin((dble(j)-f12)*dl-hpi)
  enddo
endif

!------------------------------------------------------------------
 !Process selected frames:
loop=0
do
  loop=loop+1
   !Read each frame of the data:
  read(11,rec=loop,iostat=iread) t,qq(1:ng,1:nt)
  if (iread .ne. 0) exit 

  if (iopt .eq. 2) then
     !Subtract f (Coriolis frequency):
    do i=1,nt
      qq(1:ng,i)=qq(1:ng,i)-cof(1:ng)
    enddo
  endif
    
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
      det=1.-yp**2-zp**2
      if (det .gt. zero) then
        xp=sqrt(det)
        xm=xp*clatc-zp*slatc
        zt=zp*clatc+xp*slatc
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
  write(22,rec=loop) qqop
enddo

end subroutine

!=========================================================

end program
