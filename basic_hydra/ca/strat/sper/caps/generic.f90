module generic

! Module contains generic subroutines for caps in an aperiodic geometry. 
!          *** These should not be modified ***

use constants

implicit none

contains 

!=======================================================================

subroutine average(qq,qavg)

! Computes the average value of a field qq and returns the result in qavg

implicit none

 !Passed variables:
double precision:: qq(0:ny,0:nxm1)
double precision:: qavg

 !Use trapezoidal rule in both directions:
qavg=dsumi*(f12*sum(qq(0,:)+qq(ny,:))+sum(qq(1:nym1,:)))
 !Note: dsumi = 1/(nx*ny)

return
end subroutine average

!=======================================================================

subroutine l1norm(qq,ql1)

! Computes the L1 norm of a field qq and returns the result in ql1:
!          ql1 = int_xmin^xmax{int_ymin^ymax{|qq| dxdy}}

implicit none

 !Passed variables:
double precision:: qq(0:ny,0:nxm1)
double precision:: ql1

 !Use trapezoidal rule in both directions:
ql1=garea*(f12*sum(abs(qq(0,:))+abs(qq(ny,:)))+ &
               sum(abs(qq(1:nym1,:))))
 !Note: garea is the grid box area, glx*gly

return
end subroutine l1norm

!=======================================================================

subroutine l2norm(qq,ql2)

! Computes the L2 norm of a field qq and returns the result in ql2:
!          ql2 = int_xmin^xmax{int_ymin^ymax{qq^2 dxdy}}

implicit none

 !Passed variables:
double precision:: qq(0:ny,0:nxm1)
double precision:: ql2

 !Use trapezoidal rule in both directions:
ql2=garea*(f12*sum(qq(0,:)**2+qq(ny,:)**2)+sum(qq(1:nym1,:)**2))
 !Note: garea is the grid box area, glx*gly

return
end subroutine l2norm

!=======================================================================

 !Main end module
end module
