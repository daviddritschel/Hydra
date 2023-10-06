module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral

 !Define quantities which need to be preserved between recontouring and evolution:

 !Gridded buoyancy & vorticity fields:
double precision:: bb(0:ny,0:nxm1),zs(0:ny,0:nxm1),zd(0:ny,0:nxm1)
double precision:: hh(0:nxm1)

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

 !Domain area & reference potential energy:
double precision:: domarea,eperef

contains 

!=======================================================================

subroutine interpol(fp,fe,xp,yp)

! Interpolates the field fp(iy,ix) at a given set of points
! xp(iy,ix), yp(iy,ix) (given in grid units) using 
! bi-cubic Lagrange interpolation.

! Clipping is used to restrict min & max values over each grid box.
! The interpolated field is written into the array fe.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: fp(0:ny,0:nxm1),fe(0:ny,0:nxm1)
double precision:: xp(0:ny,0:nxm1),yp(0:ny,0:nxm1)
 !Local arrays:
double precision:: phix(-1:2),phiy(-1:2)
integer:: ixper(-1:nx+2)
!---------------------------------------------------------------
 !Grid box reference index used below for x periodicity:
ixper(-1)=nxm1
do ix=0,nxm1
  ixper(ix)=ix
enddo
ixper(nx  )=0
ixper(nx+1)=1
ixper(nx+2)=2

do ix=0,nxm1
  do iy=0,ny
    ix0=int(xp(iy,ix))
     !ix0: the x grid box containing the point
    px0=xp(iy,ix)-dble(ix0)

    pxm=one+px0
    px1=one-px0
    px2=two-px0
    phix(-1)=-f16*px0*px1*px2
    phix( 0)= f12*pxm*px1*px2
    phix( 1)= f12*pxm*px0*px2
    phix( 2)=-f16*pxm*px0*px1

    iy0=min(int(yp(iy,ix)),nym1)
     !iy0: the y grid box containing the point
    py0=yp(iy,ix)-dble(iy0)

    pym=one+py0
    py1=one-py0
    py2=two-py0
    phiy(-1)=-f16*py0*py1*py2
    phiy( 0)= f12*pym*py1*py2
    phiy( 1)= f12*pym*py0*py2
    phiy( 2)=-f16*pym*py0*py1

    fe(iy,ix)=zero
    do jx=-1,2
      ixi=ixper(ix0+jx)
      do jy=-1,2
        iyi=ny-abs(ny-abs(iy0+jy))
        fe(iy,ix)=fe(iy,ix)+fp(iyi,ixi)*phix(jx)*phiy(jy)
      enddo
    enddo
     !Clip function to min/max values at corners of grid box:
    ix0p1=ixper(ix0+1)
    ix0  =ixper(ix0)
    femin=min(fp(iy0,ix0),fp(iy0+1,ix0),fp(iy0,ix0p1),fp(iy0+1,ix0p1))
    femax=max(fp(iy0,ix0),fp(iy0+1,ix0),fp(iy0,ix0p1),fp(iy0+1,ix0p1))
    fe(iy,ix)=min(femax,max(femin,fe(iy,ix)))
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine combine(qq,qc,qs,qd,qavg)

! Combine contour (qc), large scale (qs), and residual (qd) fields 
! into qq, ensuring the domain average of qq = qavg.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1),qc(0:ny,0:nxm1)
double precision:: qs(0:ny,0:nxm1),qd(0:ny,0:nxm1)
 !Define local array:
double precision:: wka(0:ny,0:nxm1)

 !Define q = F[qs-qc]+qc+qd where F is a low pass filter 
 !(see subroutine filter):
do ix=0,nxm1
  do iy=0,ny
    wka(iy,ix)=qs(iy,ix)-qc(iy,ix)
  enddo
enddo

call filter(wka,0,2)

do ix=0,nxm1
  do iy=0,ny
    qq(iy,ix)=wka(iy,ix)+qc(iy,ix)+qd(iy,ix)
  enddo
enddo

 !Restore domain average:
call restore(qq,qavg)

return
end subroutine
!=======================================================================

subroutine reset(qc,qs,qd,qavg)

! Resets the gridded fields qs & qd and ensures that <qs> = qavg

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qc(0:ny,0:nxm1),qs(0:ny,0:nxm1),qd(0:ny,0:nxm1)
 !Define local array:
double precision:: wka(0:ny,0:nxm1)

 !------------------------------------------------------------
 !Reset qs = q = F[qs-qc]+qc+qd, where F is a low pass filter:
do ix=0,nxm1
  do iy=0,ny
    wka(iy,ix)=qs(iy,ix)-qc(iy,ix)
  enddo
enddo

call filter(wka,0,2)

do ix=0,nxm1
  do iy=0,ny
    qs(iy,ix)=wka(iy,ix)+qc(iy,ix)
  enddo
enddo
 !Restore domain average:
call average(qs,qavg0)
qsadd=qavg-qavg0

 !------------------------------------------------------------
 !Reset qd = q-qc-F[q-qc]
do ix=0,nxm1
  do iy=0,ny
    qs(iy,ix)=qs(iy,ix)+qsadd+qd(iy,ix)
    qd(iy,ix)=qs(iy,ix)-qc(iy,ix)
    wka(iy,ix)=qd(iy,ix)
  enddo
enddo

call filter(wka,0,2)

do ix=0,nxm1
  do iy=0,ny
    qd(iy,ix)=qd(iy,ix)-wka(iy,ix)
  enddo
enddo
 !Ensure zero domain average qd at beginning of next time step:
call restore(qd,zero)

 !Recompute domain average (which evolves because of qd):
call average(qs,qavg)

return
end subroutine

!=======================================================================

subroutine filter(qq,isym,nrep)

! Performs nrep 1-2-1 filters in each direction to return a 
! low pass filtered version of the original array.
! If isym = 0, var is assumed to be symmetric across the boundaries, 
! while if isym = 1, var is assumed to be anti-symmetric.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1)
 !Define local array:
double precision:: wka(0:ny,0:nxm1)

 !Perform 1-2-1 average nrep times utilising intermediate work array:
if (isym .eq. 0) then
   !The function is symmetric across y boundaries:

  do j=1,nrep
    do iy=0,ny
      wka(iy,   0)=f12*qq(iy,   0)+f14*(qq(iy,nxm1)+qq(iy,1))
      wka(iy,nxm1)=f12*qq(iy,nxm1)+f14*(qq(iy,nxm2)+qq(iy,0))
    enddo
    do ix=1,nxm2
      do iy=0,ny
        wka(iy,ix)=f12*qq(iy,ix)+f14*(qq(iy,ix-1)+qq(iy,ix+1))
      enddo
    enddo

    do ix=0,nxm1
      qq(0, ix)=f12*(wka(0, ix)+wka(1,   ix))
      do iy=1,nym1
        qq(iy,ix)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
      enddo
      qq(ny,ix)=f12*(wka(ny,ix)+wka(nym1,ix))
    enddo
  enddo

else
   !Function is anti-symmetric across y boundaries:

  do j=1,nrep
    do iy=0,ny
      wka(iy,   0)=f12*qq(iy,   0)+f14*(qq(iy,nxm1)+qq(iy,1))
      wka(iy,nxm1)=f12*qq(iy,nxm1)+f14*(qq(iy,nxm2)+qq(iy,0))
    enddo
    do ix=1,nxm2
      do iy=0,ny
        wka(iy,ix)=f12*qq(iy,ix)+f14*(qq(iy,ix-1)+qq(iy,ix+1))
      enddo
    enddo

    do ix=0,nxm1
      qq(0, ix)=f12*wka(0, ix)
      do iy=1,nym1
        qq(iy,ix)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
      enddo
      qq(ny,ix)=f12*wka(ny,ix)
    enddo
  enddo

endif

return
end subroutine

!=======================================================================

subroutine average(qq,qavg)

! Computes the average value of a field qq and returns the result in qavg

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1)

 !Use trapezoidal rule in both directions:
qavg=zero
do ix=0,nxm1
  qavg=qavg+qq(0,ix)*confac(0,ix)+qq(ny,ix)*confac(ny,ix)
enddo
qavg=f12*qavg

do ix=0,nxm1
  do iy=1,nym1
    qavg=qavg+qq(iy,ix)*confac(iy,ix)
  enddo
enddo

qavg=qavg*garea/domarea

return
end subroutine

!=======================================================================

subroutine l1norm(qq,ql1)

! Computes the L1 norm of a field qq and returns the result in ql1:
!          ql1 = int_xmin^xmax{int_ymin^ymax{|qq| dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1)

 !Use trapezoidal rule in both directions:
ql1=zero
do ix=0,nxm1
  ql1=ql1+abs(qq(0,ix))*confac(0,ix)+abs(qq(ny,ix))*confac(ny,ix)
enddo
ql1=f12*ql1

do ix=0,nxm1
  do iy=1,nym1
    ql1=ql1+abs(qq(iy,ix))*confac(iy,ix)
  enddo
enddo

ql1=garea*ql1
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine l2norm(qq,ql2)

! Computes the L2 norm of a field qq and returns the result in ql2:
!          ql2 = int_xmin^xmax{int_ymin^ymax{qq^2 dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1)

 !Use trapezoidal rule in both directions:
ql2=zero
do ix=0,nxm1
  ql2=ql2+qq(0,ix)**2*confac(0,ix)+qq(ny,ix)**2*confac(ny,ix)
enddo
ql2=f12*ql2

do ix=0,nxm1
  do iy=1,nym1
    ql2=ql2+qq(iy,ix)**2*confac(iy,ix)
  enddo
enddo

ql2=garea*ql2
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine binorm(qq1,qq2,qqbi)

! Computes the integral of qq1*qq2 over the domain

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq1(0:ny,0:nxm1),qq2(0:ny,0:nxm1)

 !Use trapezoidal rule in both directions:
qqbi=zero
do ix=0,nxm1
  qqbi=qqbi+qq1(0,ix)*qq2(0,ix)*confac(0,ix)+qq1(ny,ix)*qq2(ny,ix)*confac(ny,ix)
enddo
qqbi=f12*qqbi

do ix=0,nxm1
  do iy=1,nym1
    qqbi=qqbi+qq1(iy,ix)*qq2(iy,ix)*confac(iy,ix)
  enddo
enddo

qqbi=garea*qqbi
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine kinetic(uu,vv,eke)

! Computes the kinetic energy  

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)

 !Use trapezoidal rule in both directions:
eke=zero
do ix=0,nxm1
  eke=eke+(uu( 0,ix)**2+vv( 0,ix)**2)*confac( 0,ix)**2 &
         +(uu(ny,ix)**2+vv(ny,ix)**2)*confac(ny,ix)**2
enddo
eke=f12*eke

do ix=0,nxm1
  do iy=1,nym1
    eke=eke+(uu(iy,ix)**2+vv(iy,ix)**2)*confac(iy,ix)**2
  enddo
enddo

eke=f12*garea*eke
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine potential(qq,epe)

! Computes the potential energy, the area integral of -Y*q where
! q = buoyancy

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed array:
double precision:: qq(0:ny,0:nxm1)

 !Use trapezoidal rule in both directions:
epe=zero
do ix=0,nxm1
  epe=epe+qq( 0,ix)*yori( 0,ix)*confac( 0,ix) &
         +qq(ny,ix)*yori(ny,ix)*confac(ny,ix)
enddo
epe=f12*epe

do ix=0,nxm1
  do iy=1,nym1
    epe=epe+qq(iy,ix)*yori(iy,ix)*confac(iy,ix)
  enddo
enddo

epe=-garea*epe-eperef
 !Note: garea is the grid box area, glx*gly and peref is a reference
 !      potential energy (computed in topo_caps.f90)

return
end subroutine

!=======================================================================

subroutine restore(qq,qavg)

! Restores the average of qq to the value qavg

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1)

 !Restore average (qavg):
call average(qq,qavg0)

qadd=qavg-qavg0
do ix=0,nxm1
  do iy=0,ny
    qq(iy,ix)=qq(iy,ix)+qadd
  enddo
enddo
 !Now qq has the correct average

return
end subroutine

!=======================================================================

end module
