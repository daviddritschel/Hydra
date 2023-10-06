!###########################################################################
!  The Ellipsoidal Barotropic Combined Lagrangian Advection Method (CLAM)
!  for 2D flow simulations on an ellipsoid of revolution.

!  Adapted from the spherical version in January 2015 by D.G. Dritschel
!###########################################################################

!       Uses a 4th-order Runge-Kutta scheme with an adapted time step
!       together with 4th-order compact differencing in latitude.

!       Input data file:
!       =================
!        qq_init.r8         initial gridded PV (absolute vorticity)

!       Output data files:
!       =================
!       Every approximate tsave units of time:
!        zz.r4              relative vorticity

!       Every tsim units of time (every "period"):
!        cont/synopsis.asc  number of contours, nodes, etc 
!         "   indexnnn      contour counters, etc at period "nnn" 
!         "   nodesnnn      contour nodes at period "nnn"
!        resi/pnnn          residual PV field at period "nnn"
!        grid/zz.r4         relative vorticity

!       Every contour regularisation (surgery):
!        complexity.asc  number of PV contours, nodes and time

!       Every time step:
!        monitor.asc        time, energy, angular momentum, enstrophy, etc.

!==========================================================================


!     The full algorithm consists of the following modules:
!        caps.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded fields 
!                        using PS method along with contours.
!----------------------------------------------------------------------------
program caps

use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .lt. tfin)

   !Obtain new PV contours:
  call recont

   !Advect fields until next recontouring or end:
  call evolve

enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit double precision(a-h,o-z)
implicit integer(i-n)

!--------------------------------------------------------------------
 !Initialise modules:
call init_spectral
call init_contours

 !No contours initially (the are built at the start of the run):
n=0
npt=0

!--------------------------------------------------------------------
 !Read initial gridded PV into qd for initial contouring:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(20,rec=1) t,qd
close(20)

 !Initialise gridded relative vorticity:
do i=1,nt
  do j=1,ng
    qs(j,i)=qd(j,i)-cof(j)
  enddo
enddo

 !FFT qs in longitude (where it remains as a semi-spectral array):
call forfft(ng,nt,qs,trig,factors) 

 !Set data save flag:
iref=-1

!--------------------------------------------------------------------
 !Open various diagnostic files:

 !Energy & other diagnostics:
open(18,file='monitor.asc',status='unknown')

 !Contour complexity (n, npt & t):
open(19,file='complexity.asc',status='unknown')

 !relative vorticity every tsim time units:
open(35,file='grid/zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !Relative vorticity approximately every tsave time units:
open(45,file='zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !Contour summary & residual PV files for fine-grid PV reconstruction:
open(80,file='cont/synopsis.asc',status='unknown')
open(83,file='cont/qd.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !Set dump counter for writing to correct record:
idump=0

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect fields until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine

!=======================================================================

subroutine recont

use congen

implicit none

!Obtain new PV contours:
write(*,*) 'Recontouring PV ...'
call recontour

write(*,'(a,i9,a,i10,a,f12.5)') '   n = ',n,'   npt = ',npt,'   t = ',t
write(19,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt

return 
end subroutine

!=======================================================================

subroutine finalise

implicit none

write(*,*) 'caps completed normally' 

 !Close output files (opened in subroutine initialise):
close(18)
close(19)
close(35)
close(45)
close(80)
close(83)

return 
end subroutine

 !End main program
end program
!=======================================================================
