module generic

! Module contains generic subroutines for casl in a doubly periodic geometry. 
!          *** These should not be modified ***

use constants

implicit none

contains 

!=======================================================================

subroutine average(qq,qavg)

! Computes the average value of a field qq and returns the result in qavg

implicit none

 !Passed variables:
double precision:: qq(ny,nx),qavg

 !Local variables:
integer:: ix,iy

qavg=zero
do ix=1,nx
  do iy=1,ny
    qavg=qavg+qq(iy,ix)
  enddo
enddo
qavg=qavg*dsumi

return
end subroutine

!=======================================================================

subroutine l1norm(qq,ql1)

! Computes the L1 norm of a field qq and returns the result in ql1:
!          ql1 = int_xmin^xmax{int_ymin^ymax{|qq| dxdy}}

implicit none

 !Passed variables:
double precision:: qq(ny,nx),ql1

 !Local variables:
integer:: ix,iy

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

implicit none

 !Passed variables:
double precision:: qq(ny,nx),ql2

 !Local variables:
integer:: ix,iy

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

subroutine contint(qq,nq,dq)

! Computes a contour interval for a field qq from (qq_max-qq_min)/n_q

implicit none

 !Passed variables:
double precision:: qq(ny,nx),dq
integer:: nq

 !Local variables:
double precision:: qqmin,qqmax
integer:: ix,iy

qqmin=qq(1,1)
qqmax=qq(1,1)
do ix=1,nx
  do iy=1,ny
    qqmin=min(qqmin,qq(iy,ix))
    qqmax=max(qqmax,qq(iy,ix))
  enddo
enddo

dq=(qqmax-qqmin)/dble(nq)

return
end subroutine

!=======================================================================

 !Main end module
end module
