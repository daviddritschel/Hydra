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

qavg=dsumi*sum(qq)

return
end subroutine

!=======================================================================

subroutine l1norm(qq,ql1)

! Computes the L1 norm of a field qq and returns the result in ql1:
!          ql1 = int_xmin^xmax{int_ymin^ymax{|qq| dxdy}}

implicit none

 !Passed variables:
double precision:: qq(ny,nx),ql1

ql1=garea*sum(abs(qq))
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

ql2=garea*sum(qq**2)
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

dq=(maxval(qq)-minval(qq))/dble(nq)

return
end subroutine

!=======================================================================

 !Main end module
end module
