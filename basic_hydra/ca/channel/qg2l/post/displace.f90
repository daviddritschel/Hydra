program displace
!  --------------------------------------------------------------------------
!  |  Calculates the relative displacement dd1 - dd1eq & dd2 - dd2eq        |
!  |  Writes rdd1 & rdd2.r4 (the latter only if mode 1 is not barotropic)   |
!  |  in the evolution subdirectory.                                        |
!  --------------------------------------------------------------------------

 !Import constants:
use constants

implicit none
double precision:: dd1eq(0:ny,0:nxm1), dd2eq(0:ny,0:nxm1), t
real:: dd1r4(0:ny,0:nxm1), dd2r4(0:ny,0:nxm1), tr4
integer:: loop, iread, ix, iy

!----------------------------------------------------------------
 !Read in equilibrium middle interface displacement (if present):
if (heating) then
  open(12,file='disp1eq.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) t,dd1eq
  close(12)

  if (.not. barot) then
    open(12,file='disp2eq.r8',form='unformatted', &
          access='direct',status='old',recl=2*nbytes)
    read(12,rec=1) t,dd2eq
    close(12)
  endif
else
  write(*,*) ' No equilibrium displacement fields to subtract! *** stopping ***'
  stop
endif

!---------------------------------------------------------------------
 !Open files containing interface displacements f_0*delta_j/(H_1+H_2):
open(33,file='evolution/dd1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
if (.not. barot) open(34,file='evolution/dd2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open files to contain the streamfunction in each mode:
open(43,file='evolution/rdd1.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
if (.not. barot) open(44,file='evolution/rdd2.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(33,rec=loop,iostat=iread) tr4,dd1r4
  if (iread .ne. 0) exit 
  if (.not. barot) read(34,rec=loop,iostat=iread) tr4,dd2r4

  t=dble(tr4)
  write(*,'(a,f13.5)') ' Processing t = ',t

  if (barot) then
    do ix=0,nxm1
      do iy=0,ny
        dd1r4(iy,ix)=dd1r4(iy,ix)-real(dd1eq(iy,ix))
      enddo
    enddo
    write(43,rec=loop) tr4,dd1r4
  else
    do ix=0,nxm1
      do iy=0,ny
        dd1r4(iy,ix)=dd1r4(iy,ix)-real(dd1eq(iy,ix))
        dd2r4(iy,ix)=dd2r4(iy,ix)-real(dd2eq(iy,ix))
      enddo
    enddo
    write(43,rec=loop) tr4,dd1r4
    write(44,rec=loop) tr4,dd2r4
  endif

enddo

 !Close files:
close(33)
close(43)
write(*,*)
if (.not. barot) then
  close(34)
  close(44)
  write(*,*) ' The relative interface displacements are in rdd1 & rdd2.r4.'
else
  write(*,*) ' The relative interface displacement is in rdd1.r4.'
endif

 !End main program
end program
!=======================================================================
