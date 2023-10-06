module generic

! Module contains generic subroutines for casl in a doubly periodic geometry. 
!          *** These should not be modified ***

use constants

implicit none

contains 

!=======================================================================

subroutine interpol(fp,fe,xp,yp)

! Interpolates the field fp(iy,ix) at a given set of points
! xp(iy,ix), yp(iy,ix) (given in grid units) using 
! bi-cubic Lagrange interpolation.

! The interpolated field is written into the array fe. 
! Note, the domain average value of fe is removed.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: fp(ny,nx),fe(ny,nx)
double precision:: xp(ny,nx),yp(ny,nx)
 !Local arrays:
double precision:: phix(-1:2),phiy(-1:2)
integer:: ixper(-1:nx+1),iyper(-1:ny+1)
!---------------------------------------------------------------
 !Grid box reference index used below for x periodicity:
ixper(-1)=nx
do ix=0,nx-1
  ixper(ix)=ix+1
enddo
ixper(nx)=1
ixper(nx+1)=2

iyper(-1)=ny
do iy=0,ny-1
  iyper(iy)=iy+1
enddo
iyper(ny)=1
iyper(ny+1)=2

do ix=1,nx
  do iy=1,ny
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

    iy0=int(yp(iy,ix))
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
        iyi=iyper(iy0+jy)
        fe(iy,ix)=fe(iy,ix)+fp(iyi,ixi)*phix(jx)*phiy(jy)
      enddo
    enddo
  enddo
enddo

 !Remove domain average:
call average(fe,favg)
do ix=1,nx
  do iy=1,ny
    fe(iy,ix)=fe(iy,ix)-favg
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine combine(qq,qc,qs,qd)

! Combines contour (qc), large scale (qs), and residual (qd) fields into 
! the full PV field (qq).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(ny,nx),qc(ny,nx),qs(ny,nx),qd(ny,nx)
 !Define local array:
double precision:: wka(ny,nx)

!Define q = F[qs-qc]+qc+qd
!where F is a low pass filter (see subroutine filter)

do ix=1,nx
  do iy=1,ny
    wka(iy,ix)=qs(iy,ix)-qc(iy,ix)
  enddo
enddo

call filter(wka,2)

do ix=1,nx
  do iy=1,ny
    qq(iy,ix)=wka(iy,ix)+qc(iy,ix)+qd(iy,ix)
  enddo
enddo

return
end subroutine
!=======================================================================

subroutine reset(qc,qs,qd)

! Resets the gridded fields qs & qd

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qc(ny,nx),qs(ny,nx),qd(ny,nx)
 !Define local array:
double precision:: wka(ny,nx)

 !------------------------------------------------------------
 !Reset qs = q = F[qs-qc]+qc+qd, where F is a low pass filter:
do ix=1,nx
  do iy=1,ny
    wka(iy,ix)=qs(iy,ix)-qc(iy,ix)
  enddo
enddo

call filter(wka,2)

 !Reset qd = q-qc-F[q-qc]
do ix=1,nx
  do iy=1,ny
    qs(iy,ix)=wka(iy,ix)+qc(iy,ix)+qd(iy,ix)
    qd(iy,ix)=qs(iy,ix)-qc(iy,ix)
    wka(iy,ix)=qd(iy,ix)
  enddo
enddo

call filter(wka,2)

do ix=1,nx
  do iy=1,ny
    qd(iy,ix)=qd(iy,ix)-wka(iy,ix)
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine filter(qq,nrep)

! Performs nrep 1-2-1 filters in each direction to return a 
! low pass filtered version of the original array.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(ny,nx)
 !Define local array:
double precision:: wka(ny,nx)

 !Perform 1-2-1 average nrep times utilising intermediate work array:
do j=1,nrep
  do iy=1,ny
    wka(iy, 1)=f12*qq(iy, 1)+f14*(qq(iy,  nx)+qq(iy,2))
    wka(iy,nx)=f12*qq(iy,nx)+f14*(qq(iy,nxm1)+qq(iy,1))
  enddo
  do ix=2,nxm1
    do iy=1,ny
      wka(iy,ix)=f12*qq(iy,ix)+f14*(qq(iy,ix-1)+qq(iy,ix+1))
    enddo
  enddo

  do ix=1,nx
    qq(1, ix)=f12*wka(1 ,ix)+f14*(wka(ny  ,ix)+wka(2,ix))
    qq(ny,ix)=f12*wka(ny,ix)+f14*(wka(nym1,ix)+wka(1,ix))
  enddo
  do ix=1,nx 
    do iy=2,nym1
      qq(iy,ix)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
    enddo
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine average(qq,qavg)

! Computes the average value of a field qq and returns the result in qavg

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(ny,nx)

qavg=zero
do ix=1,nx
  do iy=1,ny
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
double precision:: qq(ny,nx)

ql1=zero
do ix=1,nx
  do iy=1,ny
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
double precision:: qq(ny,nx)

ql2=zero
do ix=1,nx
  do iy=1,ny
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
double precision:: qq1(ny,nx),qq2(ny,nx)

qqbi=zero
do ix=1,nx
  do iy=1,ny
    qqbi=qqbi+qq1(iy,ix)*qq2(iy,ix)
  enddo
enddo

qqbi=garea*qqbi
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine contint(qq,nq,dq)

! Computes a contour interval for a field qq from sqrt(<q^4>/<q^2>)/n_q

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed array:
double precision:: qq(ny,nx)

ql2=zero
ql4=zero
do ix=1,nx
  do iy=1,ny
    qqsq=qq(iy,ix)**2
    ql2=ql2+qqsq
    ql4=ql4+qqsq**2
  enddo
enddo

dq=sqrt(ql4/ql2)/dble(nq)

return
end subroutine

!=======================================================================

 !Main end module
end module
