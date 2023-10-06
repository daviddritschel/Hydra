!#########################################################################
!                The Spherical Single-Layer Shallow-Water
!               Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!      Code redeveloped in June 2020 by D G Dritschel @ St Andrews.
!      Revised 11 February 2021 by dgd.
!      Equivalent Barotropic formulation added 6 January 2022 by dgd.

!      This code simulates the unforced Shallow-Water Equations (SWE)
!      or Equivalent Barotropic Equations (EBE), with R/c_p < 1,
!      in variables (q,delta,gamma), where q is the potential vorticity,
!      delta is the velocity divergence, and gamma is the acceleration 
!      divergence (called ageostrophic vorticity).  Note, we keep the
!      same definition of gamma in both the SWE and EBE.

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

   !Obtain new PV contours:
  call recont

   !Advect PV and other fields until next recontouring or end:
  call evolve

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
double precision:: qt(ng,nt)
integer:: i  

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

 !If topographic forcing is present, read in initial state:
if (forcing) then
  open(11,file='bb_init.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,bb
  close(11)
endif

!----------------------------------------------------------------------
 !Read in gridded height anomaly, h.  This helps below for defining
 !the initial fields (h is used as a guess):
open(11,file='hh_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,hh
close(11)

!----------------------------------------------------------------------
 !Read in gridded PV, qr = (zeta+f)/(1+h), where zeta is the relative
 !vorticity and h is the dimensionless height anomaly:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qr
close(11)

 !Copy into qs for proper start (see subroutine init in evolution.f90):
qs=qr

!----------------------------------------------------------------------
 !Read in gridded divergence, ds:
open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,ds
close(11)

 !Ensure domain average is zero:
call zeroavg(ds)

 !De-alias and convert to semi-spectral space as ds:
call dealias(ds)

!----------------------------------------------------------------------
 !Read in gridded acceleration divergence, gs:
open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,gs
close(11)

 !Ensure domain average is zero:
call zeroavg(gs)

 !De-alias and convert to semi-spectral space as gs:
call dealias(gs)

!----------------------------------------------------------------------
 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv),
 !relative vorticity (zz) and adjusted PV anomaly consistent with 
 !zero domain averaged zz:

 !Define PV anomaly (qt) needed for inversion below:
do i=1,nt
  qt(:,i)=qs(:,i)-cof
enddo

 !Convert qt to semi-spectral space:
call forfft(ng,nt,qt,trig,factors) 

call main_invert(qt,ds,gs,hh,uu,vv,qq,zz)
 !Note: qt, ds & gs are in semi-spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

!----------------------------------------------------------------------
 !Initially there are no contours (they are built from the gridded PV):
n=0
npt=0

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='contours/complexity.asc',status='replace')
open(15,file='evolution/ecomp.asc',status='replace')
open(16,file='evolution/ro-fr-hm.asc',status='replace')

 !Open files for coarse grid saves:
open(31,file='evolution/qq.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(32,file='evolution/dd.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(33,file='evolution/gg.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(34,file='evolution/hh.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(35,file='evolution/zz.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
if (forcing) then
   !Save topographic forcing for use in post-processing (e.g. dgbal):
  open(36,file='evolution/bb.r4',form='unformatted',access='direct', &
        status='replace',recl=nbytes)
endif

 !Open files for 1d longitudinal spectra (averaged over cos(latitude)):
open(51,file='spectra/zspec.asc',status='replace')
open(52,file='spectra/dspec.asc',status='replace')
open(53,file='spectra/gspec.asc',status='replace')
open(54,file='spectra/hspec.asc',status='replace')

 !Open files for contour writes:
open(80,file='contours/qqsynopsis.asc',status='replace')
open(83,file='contours/qqresi.r4',form='unformatted',access='direct', &
                                status='replace',recl=nbytes)

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
