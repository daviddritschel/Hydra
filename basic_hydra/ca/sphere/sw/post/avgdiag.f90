program average

! Computes the time-averaged values of the diagnostics in the file
! evolution/ro-fr-hm.asc from t_c <= t <= t_max where t_c is specified
! and t_max is the maximum time found in the data.

use parameters

implicit none

double precision:: tc, t, afac
double precision::  romax, frmax, hmin, hmax, zrms, drms
double precision:: aromax,afrmax,ahmin,ahmax,azrms,adrms
integer:: iread, na

!---------------------------------------------------------------------
! Open data file:
open(15,file='evolution/ro-fr-hm.asc',status='old')

write(*,*) 'Time from which to compute the time averages?'
read(*,*) tc
write(*,*)

! Initialise:
aromax=0.d0
afrmax=0.d0
ahmin=0.d0
ahmax=0.d0
azrms=0.d0
adrms=0.d0
na=0
  
! Read data and process:
do
  iread=0
   !Read Ro_max, Fr_max, h_min, h_max, zeta_rms, delta_rms:
  read(15,*,iostat=iread) t,romax,frmax,hmin,hmax,zrms,drms
  if (iread .ne. 0) exit 

  if (t+1.d-6 >= tc) then
    ! Accumulate time average:
    aromax=aromax+romax
    afrmax=afrmax+frmax
    ahmin=ahmin+hmin
    ahmax=ahmax+hmax
    azrms=azrms+zrms
    adrms=adrms+drms
    na=na+1
  endif

enddo

! Finalise averages:
afac=1.d0/dble(na)
aromax=afac*aromax
afrmax=afac*afrmax
ahmin=afac*ahmin
ahmax=afac*ahmax
azrms=afac*azrms
adrms=afac*adrms

! Close files:
close(51)
close(61)

write(*,'(a,f12.5)') ' The final time found is ',t
write(*,*)
write(*,*) ' The time averaged diagnotics are as follows:'
write(*,*)
write(*,*) '     Ro_max      Fr_max       h_min       h_max    zeta_rms   delta_rms'
write(*,'(6(1x,f11.6))') aromax,afrmax,ahmin,ahmax,azrms,adrms

! Record info:
open(16,file='evolution/avgdiag.asc',status='replace')
write(16,'(2(1x,f10.3),6(1x,f11.6))') tc,t,aromax,afrmax,ahmin,ahmax,azrms,adrms

! End main program
end program average
