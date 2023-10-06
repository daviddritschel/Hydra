module generic

! Module contains generic subroutines for the casl algorithm 
! in a periodic channel geometry. 

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
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine combine(qq,qc,qs,qd,qavg)

! Combines contour (qc), large scale (qs), and residual (qd) fields 
! into qq, ensuring the domain average of qq = qavg.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1),qc(0:ny,0:nxm1)
double precision:: qs(0:ny,0:nxm1),qd(0:ny,0:nxm1)
 !Define local arrays:
double precision:: wka(0:ny,0:nxm1)

!Define q = F[qs-qc]+qc+qd
!where F is a low pass filter (see subroutine filter)
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

! Resets the gridded fields qs & qd and ensures that <qs> = qavg.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qc(0:ny,0:nxm1),qs(0:ny,0:nxm1),qd(0:ny,0:nxm1)
 !Define local array:
double precision:: wka(0:ny,0:nxm1)

 !------------------------------------------------------------
 !Reset qs = q = F[qs-qc]+qc+qd, where F is a low pass filter:
wka=qs-qc
call filter(wka,0,2)
qs=wka+qc+qd

 !------------------------------------------------------------
 !Reset qd = q-qc-F[q-qc] and restore domain average for qs:
call restore(qs,qavg)
qd=qs-qc
wka=qd
 !Note: qd has a zero average by construction since the domain
 !      average of qc, like qs, is qavg.

call filter(wka,0,2)
qd=qd-wka

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

subroutine average(var,vavg)

! Computes the average value of a 2D field var and returns 
! the result in vavg

implicit none

 !Passed variables:
double precision:: var(0:ny,0:nxm1),vavg
 !Local variables:
integer:: ix,iy

 !Use trapezoidal rule in both directions:
vavg=zero
do ix=0,nxm1
  vavg=vavg+var(0,ix)+var(ny,ix)
enddo
vavg=f12*vavg

do ix=0,nxm1
  do iy=1,nym1
    vavg=vavg+var(iy,ix)
  enddo
enddo

vavg=vavg/dble(nx*ny)

return
end subroutine

!=======================================================================

subroutine restore(var,vavg)

! Restores the domain-average value of a 2D field var to vavg

implicit none

 !Passed variables:
double precision:: var(0:ny,0:nxm1),vavg
 !Local variables:
double precision:: vavg0,vadd
integer:: ix,iy

!-----------------------------------------
 !Use trapezoidal rule in both directions:
vavg0=zero
do ix=0,nxm1
  vavg0=vavg0+var(0,ix)+var(ny,ix)
enddo
vavg0=f12*vavg0

do ix=0,nxm1
  do iy=1,nym1
    vavg0=vavg0+var(iy,ix)
  enddo
enddo

vavg0=vavg0/dble(nx*ny)
vadd=vavg-vavg0
var=var+vadd

return
end subroutine

!=======================================================================

subroutine l1norm(var,vl1)

! Computes the L1 norm of a 2D field var and returns the result in vl1:
!          vl1 = int_xmin^xmax{int_ymin^ymax{|var| dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: var(0:ny,0:nxm1)

 !Use trapezoidal rule in both directions:
vl1=zero
do ix=0,nxm1
  vl1=vl1+abs(var(0,ix))+abs(var(ny,ix))
enddo
vl1=f12*vl1

do ix=0,nxm1
  do iy=1,nym1
    vl1=vl1+abs(var(iy,ix))
  enddo
enddo

vl1=garea*vl1
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine l2norm(var,vl2)

! Computes the L2 norm of a 2D field var and returns the result in vl2:
!          vl2 = int_xmin^xmax{int_ymin^ymax{var^2 dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: var(0:ny,0:nxm1)

 !Use trapezoidal rule in both directions:
vl2=zero
do ix=0,nxm1
  vl2=vl2+var(0,ix)**2+var(ny,ix)**2
enddo
vl2=f12*vl2

do ix=0,nxm1
  do iy=1,nym1
    vl2=vl2+var(iy,ix)**2
  enddo
enddo

vl2=garea*vl2
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

 !Main end module
end module
