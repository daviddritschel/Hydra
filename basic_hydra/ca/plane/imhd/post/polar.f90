program polar
!==========================================================================
!   Interpolates (A+B_0*y)_theta, j & q in polar coordinates over 
!   R-a < r < R+a and -pi < theta < pi, and also computes the 
!   projection on cos(theta) and sin(theta) (creating functions of r only).

!   Here, theta is a rotated angle, rotating at the rate omega = z0/2,
!   where z0 is the maximum vorticity in the initially circular patch.

!   R is specified below (should be consistent with the script 
!   new-ellipse), and a = C*l_eta, where l_eta = R*delta, delta =
!   sqrt(eta*z0)/u0, eta is the magnetic diffusivity, and C is an
!   O(1) constant specified below.

!   Note, (A+B_0*y)_theta is scaled by u0*rr, j & q by z0, where 
!   u0 = z0*rr/2 is the maximum velocity at the edge of a circular vortex.

!   Adapted from cs/mhdclam/sources/polar.F by D G Dritschel on 8/11/2017
!==========================================================================

 !Import constants, parameters and common arrays needed for FFTs etc:
use constants
use spectral

implicit none

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

! Flow specific parameters:
double precision,parameter:: omega=twopi,z0=two*omega,r0=5.0*pi/32.d0
double precision,parameter:: u0=r0*omega,atsca=one/(u0*r0),rjsca=one/z0

! Width of layer in diffusive lengths:
double precision,parameter:: cmult=5.d0

! Number of radial and azimuthal intervals for output:
integer,parameter:: nrr=128,nth=nx/2
double precision,parameter:: dth=twopi/dble(nth),scalef=dth/pi

! Work arrays:
double precision:: aa(ny,nx),ax(ny,nx),ay(ny,nx)
double precision:: cc(ny,nx),qq(ny,nx)
double precision:: ss(ny,nx),vtmp(ny,nx)
double precision:: at(nrr,nth),rj(nrr,nth),zz(nrr,nth)
double precision:: atcos(nrr),atsin(nrr)
double precision:: rjcos(nrr),rjsin(nrr)
double precision:: zzcos(nrr),zzsin(nrr)
double precision:: rleta,rmin,rmax,drr
double precision:: theta0,xg,yg,th,cth,sth,rr,xi,zzbar
double precision:: xx,px,pxc,yy,py,pyc
real:: qqr4(ny,nx),t
integer:: kfr,ix,iy,ith,irr,ix0,ix1,iy0,iy1
character(len=3):: cper

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Viscous length assumed in ellipse script:
rleta=glx

 !Define radial range:
rmin=r0-cmult*rleta
rmax=r0+cmult*rleta
write(*,'(2(a,f9.7))') ' Here we analyse ',rmin,' < r < ',rmax
 !Radial and azimuthal grid divisions:
drr=(rmax-rmin)/dble(nrr)

write(*,*)
write(*,*) ' Frame to analyse (use 0 for the first one)?'
read(*,*) kfr
write(cper(1:3),'(i3.3)') kfr

!---------------------------------------------------------------
 !Read data and process:
open(31,file='evolution/aa.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(31,rec=kfr+1) t,qqr4
close(31)
aa=dble(qqr4)

open(31,file='evolution/jj.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(31,rec=kfr+1) t,qqr4
close(31)
cc=dble(qqr4)

open(31,file='evolution/qq.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(31,rec=kfr+1) t,qqr4
close(31)
qq=dble(qqr4)

write(*,'(a,f12.5)') '  *** Processing t = ',t
theta0=omega*t-pi

!----------------------------------------------------------------
 !Compute (A+B_0*y)_theta.  First transform aa as ss (spectral):
call ptospc(nx,ny,aa,ss,xfactors,yfactors,xtrig,ytrig)

 !Compute x & y derivatives spectrally and transform back as ax & ay:
call xderiv(nx,ny,hrkx,ss,vtmp)
call spctop(nx,ny,vtmp,ax,xfactors,yfactors,xtrig,ytrig)
call yderiv(nx,ny,hrky,ss,vtmp)
call spctop(nx,ny,vtmp,ay,xfactors,yfactors,xtrig,ytrig)

 !Compute azimuthal derivative of (A+B_0*y) and scale by u0*r0:
do ix=1,nx
  xg=xmin+glx*dble(ix-1)
  do iy=1,ny
    yg=ymin+gly*dble(iy-1)
    aa(iy,ix)=atsca*((ay(iy,ix)+b0)*xg-ax(iy,ix)*yg)
    cc(iy,ix)=rjsca*cc(iy,ix)
    qq(iy,ix)=rjsca*qq(iy,ix)
  enddo
enddo
 !Above, aa is re-used for (A+B_0*y)_theta/(u0*r0)

 !Next interpolate (A+B_0*y)_theta/(u0*r0), j/z0 & q/z0 in polar coordinates:
do ith=1,nth
  th=theta0+dth*(dble(ith)-f12)
  cth=cos(th)
  sth=sin(th)
  do irr=1,nrr
    rr=rmin+drr*(dble(irr)-f12)

    xx=glxi*(rr*cth+pi)
    ix0=1+int(xx)
    pxc=dble(ix0)-xx
    px=one-pxc
    ix1=ix0+1

    yy=glyi*(rr*sth+pi)
    iy0=1+int(yy)
    pyc=dble(iy0)-yy
    py=one-pyc
    iy1=iy0+1

    at(irr,ith)=pyc*(pxc*aa(iy0,ix0)+px*aa(iy0,ix1)) &
                +py*(pxc*aa(iy1,ix0)+px*aa(iy1,ix1))

    rj(irr,ith)=pyc*(pxc*cc(iy0,ix0)+px*cc(iy0,ix1)) &
                +py*(pxc*cc(iy1,ix0)+px*cc(iy1,ix1))

    zz(irr,ith)=pyc*(pxc*qq(iy0,ix0)+px*qq(iy0,ix1)) &
                +py*(pxc*qq(iy1,ix0)+px*qq(iy1,ix1))
  enddo
enddo

 !Remove zonal mean value of zz:
do irr=1,nrr
  zzbar=sum(zz(irr,:))/dble(nth)
  zz(irr,:)=zz(irr,:)-zzbar
enddo

 !Write full fields:
open(41,file='evolution/at'//cper//'_polar.r4',status='replace', &
      access='direct',form='unformatted',recl=4*(1+nth*nrr))
write(41,rec=1) t,real(at)
close(41)

open(41,file='evolution/jj'//cper//'_polar.r4',status='replace', &
      access='direct',form='unformatted',recl=4*(1+nth*nrr))
write(41,rec=1) t,real(rj)
close(41)

open(41,file='evolution/zz'//cper//'_polar.r4',status='replace', &
      access='direct',form='unformatted',recl=4*(1+nth*nrr))
write(41,rec=1) t,real(zz)
close(41)

 !Now project data on cos(theta) and sin(theta):
atcos=zero
atsin=zero
rjcos=zero
rjsin=zero
zzcos=zero
zzsin=zero

do ith=1,nth
  th=dth*(dble(ith)-f12)-pi
  cth=cos(th)
  sth=sin(th)
  atcos(:)=atcos(:)+at(:,ith)*cth
  atsin(:)=atsin(:)+at(:,ith)*sth
  rjcos(:)=rjcos(:)+rj(:,ith)*cth
  rjsin(:)=rjsin(:)+rj(:,ith)*sth
  zzcos(:)=zzcos(:)+zz(:,ith)*cth
  zzsin(:)=zzsin(:)+zz(:,ith)*sth
enddo

atcos=scalef*atcos
atsin=scalef*atsin
rjcos=scalef*rjcos
rjsin=scalef*rjsin
zzcos=scalef*zzcos
zzsin=scalef*zzsin

 !Write cross sections:
open(35,file='evolution/at-jj-zz'//cper//'_cross.dat',status='replace')
do irr=1,nrr
  rr=drr*(dble(irr)-f12)+rmin
   !radial coordinate scaled on the viscous length:
  xi=(rr-r0)/rleta
  write(35,'(1x,f12.7,6(1x,e14.7))') xi,atcos(irr),atsin(irr), &
                                        rjcos(irr),rjsin(irr), &
                                        zzcos(irr),zzsin(irr)
enddo
close(35)

write(*,*)
write(*,*) ' All done.'
write(*,*)
write(*,*) ' Image (A+B_0*y)_theta/(u_0*r0) using'
write(*,'(a,i4,a,i3,a)') 'dataview -ndim ',nth,' ',nrr,' evolution/at'//cper//'_polar.r4'
write(*,*)
write(*,*) ' Image j/z0 using'
write(*,'(a,i4,a,i3,a)') 'dataview -ndim ',nth,' ',nrr,' evolution/jj'//cper//'_polar.r4'
write(*,*)
write(*,*) ' Image q/z0 using'
write(*,'(a,i4,a,i3,a)') 'dataview -ndim ',nth,' ',nrr,' evolution/zz'//cper//'_polar.r4'

write(*,*)
write(*,*) ' For the projections on cos(theta) and sin(theta),'
write(*,*) ' see the file evolution/at-jj-zz'//cper//'_cross.dat'

return
end subroutine

 !End main program
end program
!=======================================================================
