program extract
!=========================================================================!
!  Extracts the gridded buoyancy field at a chosen time and writes the
!  data to a file.
!=========================================================================!

 !Import constants and parameters:
use constants

implicit none

real:: qqr4(ny,nx)
real:: t
integer:: loop

!--------------------------------------------------------------
write(*,*) ' Enter the time to extract:'
read(*,*) t
loop=nint(t/tgsave)+1

open(31,file='bb.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(31,rec=loop) t,qqr4
close(31)

open(41,file='qq_extract.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(41,rec=1) zero,dble(qqr4)
close(41)

write(*,'(a,f7.2,a)') ' Buoyancy at t = ',t,' written to qq_extract.r8'

end program extract
