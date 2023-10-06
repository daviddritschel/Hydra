module generic

! Module contains generic subroutines for casl in an aperiodic geometry. 
!          *** These should not be modified ***

use constants

implicit none

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
double precision:: fp(0:ny,0:nx),xp(0:ny,0:nx),yp(0:ny,0:nx),fe(0:ny,0:nx)

 !Local arrays:
double precision:: phix(-1:2),phiy(-1:2)

!---------------------------------------------------------------
do ix=0,nx
  do iy=0,ny
    ix0=min(int(xp(iy,ix)),nxm1)
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
      ixi=nx-abs(nx-abs(ix0+jx))
      do jy=-1,2
        iyi=ny-abs(ny-abs(iy0+jy))
        fe(iy,ix)=fe(iy,ix)+fp(iyi,ixi)*phix(jx)*phiy(jy)
      enddo
    enddo
     !Clip function to min/max values at corners of grid box:
    femin=min(fp(iy0,ix0),fp(iy0+1,ix0),fp(iy0,ix0+1),fp(iy0+1,ix0+1))
    femax=max(fp(iy0,ix0),fp(iy0+1,ix0),fp(iy0,ix0+1),fp(iy0+1,ix0+1))
    fe(iy,ix)=min(femax,max(femin,fe(iy,ix)))
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine interpolqe(fp,fe,xp,yp)

! Interpolates the field fp(iy,ix) at a given set of points
! xp(iy,ix), yp(iy,ix) (given in grid units) using 
! bi-cubic Lagrange interpolation except near the boundaries, 
! where bi-quadratic interpolation is used.

! Clipping is used to restrict min & max values over each grid box.
! The interpolated field is written into the array fe.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: fp(0:ny,0:nx),fe(0:ny,0:nx),xp(0:ny,0:nx),yp(0:ny,0:nx)

 !Local arrays:
double precision:: femin(0:ny,0:nx),femax(0:ny,0:nx)
double precision:: phix(-1:2),phiy(-1:2)

!---------------------------------------------------------------
 !Define min & max values over each grid box for clipping below:
do ix=0,nxm1
  do iy=0,nym1
    femin(iy,ix)=min(fp(iy,ix),fp(iy+1,ix),fp(iy,ix+1),fp(iy+1,ix+1))
    femax(iy,ix)=max(fp(iy,ix),fp(iy+1,ix),fp(iy,ix+1),fp(iy+1,ix+1))
  enddo
enddo

do ix=0,nx
  do iy=0,ny
    ix0=min(int(xp(iy,ix)),nxm1)
     !ix0: the x grid box containing the point
    px0=xp(iy,ix)-dble(ix0)
    if (ix0 .eq. 0) then
       !Left edge; use quadratic interpolation:
      jxbeg=0
      jxend=2
      px1=one-px0
      px2=two-px0
      phix( 0)= f12*px1*px2
      phix( 1)=     px0*px2
      phix( 2)=-f12*px0*px1
    else if (ix0 .eq. nxm1) then
       !Right edge; use quadratic interpolation:
      jxbeg=-1
      jxend=1
      pxm=one+px0
      px1=one-px0
      phix(-1)=-f12*px0*px1
      phix( 0)=     pxm*px1
      phix( 1)= f12*pxm*px0
    else
       !Interior; use cubic interpolation:
      jxbeg=-1
      jxend=2
      pxm=one+px0
      px1=one-px0
      px2=two-px0
      phix(-1)=-f16*px0*px1*px2
      phix( 0)= f12*pxm*px1*px2
      phix( 1)= f12*pxm*px0*px2
      phix( 2)=-f16*pxm*px0*px1
    endif

    iy0=min(int(yp(iy,ix)),nym1)
     !iy0: the y grid box containing the point
    py0=yp(iy,ix)-dble(iy0)

    if (iy0 .eq. 0) then
       !Bottom edge; use quadratic interpolation:
      jybeg=0
      jyend=2
      py1=one-py0
      py2=two-py0
      phiy( 0)= f12*py1*py2
      phiy( 1)=     py0*py2
      phiy( 2)=-f12*py0*py1
    else if (iy0 .eq. nym1) then
       !Top edge; use quadratic interpolation:
      jybeg=-1
      jyend=1
      pym=one+py0
      py1=one-py0
      phiy(-1)=-f12*py0*py1
      phiy( 0)=     pym*py1
      phiy( 1)= f12*pym*py0
    else
       !Interior; use cubic interpolation:
      jybeg=-1
      jyend=2
      pym=one+py0
      py1=one-py0
      py2=two-py0
      phiy(-1)=-f16*py0*py1*py2
      phiy( 0)= f12*pym*py1*py2
      phiy( 1)= f12*pym*py0*py2
      phiy( 2)=-f16*pym*py0*py1
    endif

    fe(iy,ix)=zero
    do jx=jxbeg,jxend
      do jy=jybeg,jyend
        fe(iy,ix)=fe(iy,ix)+fp(iy0+jy,ix0+jx)*phix(jx)*phiy(jy)
      enddo
    enddo
     !Clip function to min/max values at corners of grid box:
    fe(iy,ix)=min(femax(iy0,ix0),max(femin(iy0,ix0),fe(iy,ix)))
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
double precision:: qq(0:ny,0:nx),qc(0:ny,0:nx),qs(0:ny,0:nx),qd(0:ny,0:nx)
 !Define local array:
double precision:: wka(0:ny,0:nx)

!Define q = F[qs-qc]+qc+qd
!where F is a low pass filter (see subroutine filter)

do ix=0,nx
  do iy=0,ny
    wka(iy,ix)=qs(iy,ix)-qc(iy,ix)
  enddo
enddo

call filter(wka,0,2)

do ix=0,nx
  do iy=0,ny
    qq(iy,ix)=wka(iy,ix)+qc(iy,ix)
  enddo
enddo

 !Restore domain average:
call average(qq,qavg0)
qqadd=qavg-qavg0
do ix=0,nx
  do iy=0,ny
    qq(iy,ix)=qq(iy,ix)+qqadd+qd(iy,ix)
  enddo
enddo

return
end subroutine
!=======================================================================

subroutine reset(qc,qs,qd,qavg)

! Resets the gridded fields qs & qd and ensures that <qs> = qavg

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qc(0:ny,0:nx),qs(0:ny,0:nx),qd(0:ny,0:nx)
 !Define local array:
double precision:: wka(0:ny,0:nx)

 !------------------------------------------------------------
 !Reset qs = q = F[qs-qc]+qc+qd, where F is a low pass filter:
do ix=0,nx
  do iy=0,ny
    wka(iy,ix)=qs(iy,ix)-qc(iy,ix)
  enddo
enddo

call filter(wka,0,2)

do ix=0,nx
  do iy=0,ny
    qs(iy,ix)=wka(iy,ix)+qc(iy,ix)
  enddo
enddo
 !Restore domain average:
call average(qs,qavg0)
qsadd=qavg-qavg0

 !------------------------------------------------------------
 !Reset qd = q-qc-F[q-qc]
do ix=0,nx
  do iy=0,ny
    qs(iy,ix)=qs(iy,ix)+qsadd+qd(iy,ix)
    qd(iy,ix)=qs(iy,ix)-qc(iy,ix)
    wka(iy,ix)=qd(iy,ix)
  enddo
enddo

call filter(wka,0,2)

do ix=0,nx
  do iy=0,ny
    qd(iy,ix)=qd(iy,ix)-wka(iy,ix)
  enddo
enddo

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
double precision:: qq(0:ny,0:nx)
 !Define local array:
double precision:: wka(0:ny,0:nx)

 !Perform 1-2-1 average nrep times utilising intermediate work array:
if (isym .eq. 0) then
   !The function is symmetric across boundaries:

  do j=1,nrep
    do iy=0,ny
      wka(iy, 0)=f12*(qq(iy, 0)+qq(iy,   1))
      wka(iy,nx)=f12*(qq(iy,nx)+qq(iy,nxm1))
    enddo
    do ix=1,nxm1
      do iy=0,ny
        wka(iy,ix)=f12*qq(iy,ix)+f14*(qq(iy,ix-1)+qq(iy,ix+1))
      enddo
    enddo

    do ix=0,nx
      qq(0, ix)=f12*(wka(0, ix)+wka(1,   ix))
      do iy=1,nym1
        qq(iy,ix)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
      enddo
      qq(ny,ix)=f12*(wka(ny,ix)+wka(nym1,ix))
    enddo
  enddo

else
   !Function is anti-symmetric across boundaries:

  do j=1,nrep
    do iy=0,ny
      wka(iy, 0)=f12*qq(iy, 0)
      wka(iy,nx)=f12*qq(iy,nx)
    enddo
    do ix=1,nxm1
      do iy=0,ny
        wka(iy,ix)=f12*qq(iy,ix)+f14*(qq(iy,ix-1)+qq(iy,ix+1))
      enddo
    enddo

    do ix=0,nx
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
double precision:: qq(0:ny,0:nx)

 !Use trapezoidal rule in both directions:
qavg=zero
do ix=1,nxm1
  qavg=qavg+qq(0,ix)+qq(ny,ix)
enddo
do iy=1,nym1
  qavg=qavg+qq(iy,0)+qq(iy,nx)
enddo
qavg=f12*qavg+f14*(qq(0,0)+qq(ny,0)+qq(0,nx)+qq(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    qavg=qavg+qq(iy,ix)
  enddo
enddo

qavg=qavg/dble(nx*ny)

return
end subroutine

!=======================================================================

subroutine l1norm(qq,ql1)

! Computes the L1 norm of a field qq and returns the result in ql1:
!          ql1 = int_xmin^xmax{int_ymin^ymax{|qq| dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nx)

 !Use trapezoidal rule in both directions:
ql1=zero
do ix=1,nxm1
  ql1=ql1+abs(qq(0,ix))+abs(qq(ny,ix))
enddo
do iy=1,nym1
  ql1=ql1+abs(qq(iy,0))+abs(qq(iy,nx))
enddo
ql1=f12*ql1+f14*(abs(qq(0,0))+abs(qq(ny,0))+abs(qq(0,nx))+abs(qq(ny,nx)))

do ix=1,nxm1
  do iy=1,nym1
    ql1=ql1+abs(qq(iy,ix))
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
double precision:: qq(0:ny,0:nx)

 !Use trapezoidal rule in both directions:
ql2=zero
do ix=1,nxm1
  ql2=ql2+qq(0,ix)**2+qq(ny,ix)**2
enddo
do iy=1,nym1
  ql2=ql2+qq(iy,0)**2+qq(iy,nx)**2
enddo
ql2=f12*ql2+f14*(qq(0,0)**2+qq(ny,0)**2+qq(0,nx)**2+qq(ny,nx)**2)

do ix=1,nxm1
  do iy=1,nym1
    ql2=ql2+qq(iy,ix)**2
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
double precision:: qq1(0:ny,0:nx),qq2(0:ny,0:nx)

 !Use trapezoidal rule in both directions:
qqbi=zero
do ix=1,nxm1
  qqbi=qqbi+qq1(0,ix)*qq2(0,ix)+qq1(ny,ix)*qq2(ny,ix)
enddo
do iy=1,nym1
  qqbi=qqbi+qq1(iy,0)*qq2(iy,0)+qq1(iy,nx)*qq2(iy,nx)
enddo
qqbi=f12*qqbi+f14*(qq1(0,0) *qq2(0,0) +qq1(ny,0) *qq2(ny,0) &
                & +qq1(0,nx)*qq2(0,nx)+qq1(ny,nx)*qq2(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    qqbi=qqbi+qq1(iy,ix)*qq2(iy,ix)
  enddo
enddo

qqbi=garea*qqbi
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

 !Main end module
end module
