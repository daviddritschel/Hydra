program slope
!  --------------------------------------------------------------------------
!  |  Calculates the domain average slope of the interface displacements    |
!  |  Writes mean-slope.asc in the diagnostics subdirectory.                |
!  --------------------------------------------------------------------------

 !Import constants:
use constants

implicit none
real:: dd1(0:ny,0:nxm1), dd2(0:ny,0:nxm1), t
real:: sl1, sl2
integer:: loop, iread, ix

!---------------------------------------------------------------------
 !Open files containing interface displacements f_0*delta_j/(H_1+H_2):
open(33,file='evolution/dd1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
if (.not. barot) open(34,file='evolution/dd2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open files to contain the mean slope:
open(55,file='diagnostics/mean-slope.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(33,rec=loop,iostat=iread) t,dd1
  if (iread .ne. 0) exit 
  if (.not. barot) read(34,rec=loop,iostat=iread) t,dd2

  write(*,'(a,f13.5)') ' Processing t = ',t

  sl1=zero
  sl2=zero

  if (barot) then
    do ix=0,nxm1
      sl1=sl1+dd1(ny,ix)-dd1(0,ix)
    enddo
    sl1=sl1/(float(nx)*elly)
  else
    do ix=0,nxm1
      sl1=sl1+dd1(ny,ix)-dd1(0,ix)
      sl2=sl2+dd2(ny,ix)-dd2(0,ix)
    enddo
    sl1=sl1/(float(nx)*elly)
    sl2=sl2/(float(nx)*elly)
  endif

  write(55,'(f9.2,2(1x,f12.7))') t,sl1,sl2

enddo

 !Close files:
close(33)
if (.not. barot) close(34)
close(55)
  
write(*,*) ' The mean slopes vs time are in mean-slope.asc'

 !End main program
end program
!=======================================================================
