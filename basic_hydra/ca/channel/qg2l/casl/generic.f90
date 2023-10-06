module generic

! Module contains generic subroutines for the casl algorithm 
! in a multi-layer aperiodic geometry. 

!          *** These should not be modified ***

use constants

implicit none

contains 

!=======================================================================

subroutine interpol(fp,fe,xp,yp)

! Interpolates the field fp(iy,ix,iz) at a given set of points
! xp(iy,ix,iz), yp(iy,ix,iz) (given in grid units) using 
! bi-cubic Lagrange interpolation.

! The interpolated field is written into the array fe.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: fp(0:ny,0:nxm1,nz),fe(0:ny,0:nxm1,nz)
double precision:: xp(0:ny,0:nxm1,nz),yp(0:ny,0:nxm1,nz)
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

do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      ix0=int(xp(iy,ix,iz))
       !ix0: the x grid box containing the point
      px0=xp(iy,ix,iz)-dble(ix0)

      pxm=one+px0
      px1=one-px0
      px2=two-px0
      phix(-1)=-f16*px0*px1*px2
      phix( 0)= f12*pxm*px1*px2
      phix( 1)= f12*pxm*px0*px2
      phix( 2)=-f16*pxm*px0*px1

      iy0=min(int(yp(iy,ix,iz)),nym1)
       !iy0: the y grid box containing the point
      py0=yp(iy,ix,iz)-dble(iy0)

      pym=one+py0
      py1=one-py0
      py2=two-py0
      phiy(-1)=-f16*py0*py1*py2
      phiy( 0)= f12*pym*py1*py2
      phiy( 1)= f12*pym*py0*py2
      phiy( 2)=-f16*pym*py0*py1

      fe(iy,ix,iz)=zero
      do jx=-1,2
        ixi=ixper(ix0+jx)
        do jy=-1,2
          iyi=ny-abs(ny-abs(iy0+jy))
          fe(iy,ix,iz)=fe(iy,ix,iz)+fp(iyi,ixi,iz)*phix(jx)*phiy(jy)
        enddo
      enddo
    enddo
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine combine(qq,qc,qs,qd,qavg)

! Combines contour (qc), large scale (qs), and residual (qd) fields 
! into qq, ensuring the domain average of qq = qavg in each layer.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1,nz),qc(0:ny,0:nxm1,nz)
double precision:: qs(0:ny,0:nxm1,nz),qd(0:ny,0:nxm1,nz)
double precision:: qavg(nz)
 !Define local arrays:
double precision:: wka(0:ny,0:nxm1,nz)

!Define q = F[qs-qc]+qc+qd
!where F is a low pass filter (see subroutine filter)
do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      wka(iy,ix,iz)=qs(iy,ix,iz)-qc(iy,ix,iz)
    enddo
  enddo
enddo

call filter(wka,0,2)

do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      qq(iy,ix,iz)=wka(iy,ix,iz)+qc(iy,ix,iz)+qd(iy,ix,iz)
    enddo
  enddo
enddo

 !Restore domain average:
do iz=1,nz
  call average(qq(0,0,iz),qavg0)
  qqadd=qavg(iz)-qavg0
  do ix=0,nxm1
    do iy=0,ny
      qq(iy,ix,iz)=qq(iy,ix,iz)+qqadd
    enddo
  enddo
enddo

return
end subroutine
!=======================================================================

subroutine reset(qc,qs,qd,qavg)

! Resets the gridded fields qs & qd and ensures that <qs> = qavg
! in each layer

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qc(0:ny,0:nxm1,nz),qs(0:ny,0:nxm1,nz),qd(0:ny,0:nxm1,nz)
double precision:: qavg(nz)
 !Define local array:
double precision:: wka(0:ny,0:nxm1,nz)

 !------------------------------------------------------------
 !Reset qs = q = F[qs-qc]+qc+qd, where F is a low pass filter:
do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      wka(iy,ix,iz)=qs(iy,ix,iz)-qc(iy,ix,iz)
    enddo
  enddo
enddo

call filter(wka,0,2)

do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      qs(iy,ix,iz)=wka(iy,ix,iz)+qc(iy,ix,iz)+qd(iy,ix,iz)
    enddo
  enddo
enddo

 !------------------------------------------------------------
 !Reset qd = q-qc-F[q-qc] and restore domain average for qs:
do iz=1,nz
  call average(qs(0,0,iz),qavg0)
  qsadd=qavg(iz)-qavg0
  do ix=0,nxm1
    do iy=0,ny
      qs(iy,ix,iz)=qs(iy,ix,iz)+qsadd
      qd(iy,ix,iz)=qs(iy,ix,iz)-qc(iy,ix,iz)
      wka(iy,ix,iz)=qd(iy,ix,iz)
    enddo
  enddo
enddo
 !Note: qd has a zero average by construction since the domain
 !      average of qc, like qs, is qavg.

call filter(wka,0,2)

do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      qd(iy,ix,iz)=qd(iy,ix,iz)-wka(iy,ix,iz)
    enddo
  enddo
enddo

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
double precision:: qq(0:ny,0:nxm1,nz)
 !Define local array:
double precision:: wka(0:ny,0:nxm1)

 !Perform 1-2-1 average nrep times utilising intermediate work array:
if (isym .eq. 0) then
   !The function is symmetric across y boundaries:

  do j=1,nrep
    do iz=1,nz
      do iy=0,ny
        wka(iy,   0)=f12*qq(iy,   0,iz)+f14*(qq(iy,nxm1,iz)+qq(iy,1,iz))
        wka(iy,nxm1)=f12*qq(iy,nxm1,iz)+f14*(qq(iy,nxm2,iz)+qq(iy,0,iz))
      enddo
      do ix=1,nxm2
        do iy=0,ny
          wka(iy,ix)=f12*qq(iy,ix,iz)+f14*(qq(iy,ix-1,iz)+qq(iy,ix+1,iz))
        enddo
      enddo

      do ix=0,nxm1
        qq(0, ix,iz)=f12*(wka(0, ix)+wka(1,   ix))
        do iy=1,nym1
          qq(iy,ix,iz)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
        enddo
        qq(ny,ix,iz)=f12*(wka(ny,ix)+wka(nym1,ix))
      enddo
    enddo
  enddo

else
   !Function is anti-symmetric across y boundaries:

  do j=1,nrep
    do iz=1,nz
      do iy=0,ny
        wka(iy,   0)=f12*qq(iy,   0,iz)+f14*(qq(iy,nxm1,iz)+qq(iy,1,iz))
        wka(iy,nxm1)=f12*qq(iy,nxm1,iz)+f14*(qq(iy,nxm2,iz)+qq(iy,0,iz))
      enddo
      do ix=1,nxm2
        do iy=0,ny
          wka(iy,ix)=f12*qq(iy,ix,iz)+f14*(qq(iy,ix-1,iz)+qq(iy,ix+1,iz))
        enddo
      enddo

      do ix=0,nxm1
        qq(0, ix,iz)=f12*wka(0, ix)
        do iy=1,nym1
          qq(iy,ix,iz)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
        enddo
        qq(ny,ix,iz)=f12*wka(ny,ix)
      enddo
    enddo
  enddo

endif

return
end subroutine

!=======================================================================

subroutine average(var,vavg)

! Computes the average value of a 2D field var and returns 
! the result in vavg

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed array:
double precision:: var(0:ny,0:nxm1)

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

subroutine contint(qq,nqjumps,dq,qqmin,qqmax)

! Computes a contour interval for a field qq from (qq_max-qq_min)/n_qjumps

implicit none

 !Passed variables:
double precision:: qq(0:ny,0:nxm1,nz)
double precision:: dq(nz),qqmin(nz),qqmax(nz)
integer:: nqjumps

 !Local variables:
integer:: ix,iy,iz

do iz=1,nz
  qqmax(iz)=qq(0,0,iz)
  qqmin(iz)=qq(0,0,iz)
  do ix=0,nxm1
    do iy=0,ny
      qqmax(iz)=max(qqmax(iz),qq(iy,ix,iz))
      qqmin(iz)=min(qqmin(iz),qq(iy,ix,iz))
    enddo
  enddo
  dq(iz)=(qqmax(iz)-qqmin(iz))/dble(nqjumps)
enddo

return
end subroutine

!=======================================================================

double precision function rand(i)
 !Returns a random number: i is any integer

implicit double precision(a-h,o-z)
implicit integer(i-n)

call random_number(r)
rand=r

return
end function

!=======================================================================

 !Main end module
end module
