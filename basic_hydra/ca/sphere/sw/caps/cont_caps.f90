!#########################################################################
!                The Spherical Single-Layer Shallow-Water
!               Combined Lagrangian Advection Method (CLAM)

!      ==>  Continues from last time in contours/qqsynopsis.asc  <==
!#########################################################################

!      Code redeveloped in June 2020 by D G Dritschel @ St Andrews.
!      Revised 11 February 2021 by D G Dritschel @ St Andrews.
!      Continuation code adapted from caps.f90 on 18 April 2021 by DGD.

!      This code simulates the unforced Shallow-Water Equations (SWE) 
!      in variables (q,delta,gamma), where q is the potential vorticity,
!      delta is the velocity divergence, and gamma is the acceleration 
!      divergence (called ageostrophic vorticity).

!      Contour advection and generation are done internally now.  For
!      details of the method, see Dritschel & Fontane, J. Comput. Phys.
!      229, pp. 5408--5417 (2010).

!      The full algorithm consists of the following modules:

!      caps.f90      : This source - main program loop, repeats successive 
!                      calls to evolve fields and recontour;
!      parameters.f90: User defined parameters for a simulation;
!      constants.f90 : Fixed constants used throughout the other modules;
!      common.f90    : Common data preserved throughout simulation 
!                      (through recontouring--evolution cycle);
!      spectral.f90  : Fourier transform common storage and routines;
!      contours.f90  : Contour advection common storage and routines;
!      congen.f90    : Source code for contour-to-grid conversion;
!      evolution.f90 : Main time evolution module - advects gridded 
!                      fields using a PS method along with contours.
!----------------------------------------------------------------------------
program caps

use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Advect PV and other fields until next recontouring or end:
  call evolve

   !Obtain new PV contours:
  call recont

enddo

 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

call finalise

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

! Local variables:
double precision:: qt(ng,nt),td,tf
real:: qqr4(ng,nt),tr4
integer:: irec,m,i,j
character(len=3):: pind

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!----------------------------------------------------------------------
 !Read last time contours were written:
open(80,file='contours/qqsynopsis.asc',status='old',access='append')
backspace(80)
read(80,*) n,npt,t
backspace(80) !The line above will be rewritten in evolution.f90

write(*,'(a,f12.5)') ' Re-starting from t = ',t
write(*,*)

 !Read contours at this time:
irec=nint(t/tcsave)+1
write(pind(1:3),'(i3.3)') irec

open(81,file='contours/qqindex'//pind,form='unformatted',status='old')
read(81) np(1:n),i1(1:n),ind(1:n)
close(81)

open(82,file='contours/qqnodes'//pind,form='unformatted',status='old')
read(82) x(1:npt),y(1:npt),z(1:npt)
close(82)

open(83,file='contours/qqresi.r4',form='unformatted',access='direct', &
                                status='old',recl=nbytes)

 !Define and next(i) and i2(j):
do i=1,npt
  next(i)=i+1
enddo

do j=1,n
  i2(j)=i1(j)+np(j)-1
  next(i2(j))=i1(j)
enddo

!----------------------------------------------------------------------
 !Open other files and position at correct time:

open(16,file='evolution/ro-fr-hm.asc',status='old',access='append')
backspace(16)
read(16,*) td

if (abs(t-td) .lt. dt) then
   !Starting from a completed job; open files in append mode:
  open(14,file='contours/complexity.asc',status='old',access='append')
  open(15,file='evolution/ecomp.asc',status='old',access='append')
  open(51,file='spectra/zspec.asc',status='old',access='append')
  open(52,file='spectra/dspec.asc',status='old',access='append')
  open(53,file='spectra/gspec.asc',status='old',access='append')
  open(54,file='spectra/hspec.asc',status='old',access='append')

else
   !This is a restart from data which was not synced:
  tf=t-0.1d0*dt
  open(14,file='contours/complexity.asc',status='old')
  td=zero
  do while (td .lt. tf)
    read(14,*) td
  enddo
  backspace(14)

  open(15,file='evolution/ecomp.asc',status='old')
  td=zero
  do while (td .lt. tf)
    read(15,*) td
  enddo
  backspace(15)

  close(16)
  open(16,file='evolution/ro-fr-hm.asc',status='old')
  td=zero
  do while (td .lt. tf)
    read(16,*) td
  enddo
  backspace(16)

   !Open files for 1d longitudinal spectra (averaged over cos(latitude)):
  open(51,file='spectra/zspec.asc',status='old')
  open(52,file='spectra/dspec.asc',status='old')
  open(53,file='spectra/gspec.asc',status='old')
  open(54,file='spectra/hspec.asc',status='old')
  td=zero
  do
    read(51,*) td
    if (tf .gt. tf-tgsave) exit 
    read(52,*) td
    read(53,*) td
    read(54,*) td
    do m=1,ng
      read(51,*)
      read(52,*)
      read(53,*)
      read(54,*)
    enddo
  enddo
  backspace(51)
endif

 !Open files for coarse grid saves:
open(31,file='evolution/qq.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(32,file='evolution/dd.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(33,file='evolution/gg.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(34,file='evolution/hh.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(35,file='evolution/zz.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
 !Read current state:
irec=nint(t/tgsave)+1
read(31,rec=irec) tr4,qqr4
qq=dble(qqr4)
read(32,rec=irec) tr4,qqr4
ds=dble(qqr4)
read(33,rec=irec) tr4,qqr4
gs=dble(qqr4)
read(34,rec=irec) tr4,qqr4
hh=dble(qqr4)
read(35,rec=irec) tr4,qqr4
zz=dble(qqr4)

if (forcing) then
   !Topographic forcing:
  open(36,file='evolution/bb.r4',form='unformatted',access='direct', &
                               status='old',recl=nbytes)
   !Read current state:
  read(36,rec=irec) tr4,qqr4
  bb=dble(qqr4)
endif

 !Define quantities needed by evolve after leaving this subroutine:
qs=qq
 !Convert ds to semi-spectral space:
call forfft(ng,nt,ds,trig,factors) 
 !Convert gs to semi-spectral space:
call forfft(ng,nt,gs,trig,factors) 

!Define PV anomaly (qt) needed for inversion below:
do i=1,nt
  qt(:,i)=qq(:,i)-cof
enddo

 !Convert qt to semi-spectral space:
call forfft(ng,nt,qt,trig,factors) 

call main_invert(qt,ds,gs,hh,uu,vv,qq,zz)
 !Note: qt, ds & gs are in semi-spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

!----------------------------------------------------------------------
 !Define number of time steps between grid and contour saves:
ngsave=nint(tgsave/dt)
ncsave=nint(tcsave/dt)
 !*** WARNING: tgsave and tcsave should be an integer multiple of dt

return
end subroutine initialise

!=======================================================================

subroutine evolve

use evolution

implicit none

 !Advect PV until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return
end subroutine evolve

!=======================================================================

subroutine recont

use congen

implicit none

 !Obtain new PV contours:
write(*,*) 'Recontouring PV ...'
call recontour
write(*,'(a,i8,a,i9)') '   n = ',n,'   npt = ',npt

return 
end subroutine recont

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(14)
close(15)
close(16)
close(31)
close(32)
close(33)
close(34)
close(35)
if (forcing) close(36)
close(51)
close(52)
close(53)
close(54)
close(80)
close(83)

return
end subroutine finalise

 !End main program
end program caps
!=======================================================================
