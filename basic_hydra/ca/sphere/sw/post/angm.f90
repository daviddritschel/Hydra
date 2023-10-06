program average

! Computes the time average angular momentum over a specified period

use parameters

implicit none

double precision:: tc, t, dum1, dum2, dum3, angm, avangm
integer:: iread, na

!---------------------------------------------------------------------
! Open data file:
open(15,file='evolution/ecomp.asc',status='old')

write(*,*) 'Time from which to compute time-average angular momentum?'
read(*,*) tc
write(*,*)

! Initialise:
avangm=0.d0
na=0
  
! Read data and process:
do
  iread=0
  read(15,*,iostat=iread) t,dum1,dum2,dum3,angm
  if (iread .ne. 0) exit 

  if (t+1.d-6 >= tc) then
    ! Accumulate time average:
    avangm=avangm+angm
    na=na+1
  endif

enddo

! Finalise average:
avangm=avangm/dble(na)

! Close files:
close(51)
close(61)

write(*,*)
write(*,'(a,f9.5)') ' The time-average angular momentum is ',avangm

! End main program
end program average
