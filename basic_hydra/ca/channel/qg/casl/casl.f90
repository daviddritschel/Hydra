!#########################################################################
!        The Combined Lagrangian Advection Method for single-layer
!           quasi-geostrophic flow in periodic channel geometry
!#########################################################################

!        Code adapted from two-layer version by David Dritschel
!        @ St Andrews on 13 April 2023.

!        Uses contour advection and a semi-Lagrangian method for
!        the gridded PV fields.  Inversion is carried out using a
!        semi-spectral/4th-order difference method.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!        This code solves Dq/Dt = 0 where the PV is given by
!               q = beta*y + Lap(psi) - k_d^2 * psi
!        in the periodic channel -L_x/2 < x < L_x/2 ; -L_y/2 < y < L_y/2
!        (with free slip boundary conditions in y, periodic in x),
!        where psi is the streamfunction from which u = -dpsi/dy and
!        v = dpsi/dx. At the y boundaries, the mean x velocity takes
!        prescribed values (see parameters.f90 and spectral.f90).
!        Note: an additional mean x velocity may be included.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!        The full algorithm consists of the following modules:

!        casl.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        evolution.f90 : Main time evolution module;
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        congen.f90    : Source code for contour-to-grid conversion;
!        generic.f90   : Generic service routines for CASL.
!----------------------------------------------------------------------------
program casl

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

   !Advect PV until next recontouring or end:
  call evolve

enddo
 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Close all files and end: 
call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, opens and reads
! input files, then opens output files.

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: qqmin,qqmax

!--------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Initialise constants and arrays for contour advection:
call init_contours

!-----------------------------------------------------------------
 !Offsets used to determine trajectory location in the grid during
 !semi-Lagrangian advection:
do ix=0,nxm1
  xig(ix)=dble(ix)
enddo

do iy=0,ny
  yig(iy)=dble(iy)
enddo

!-----------------------------------------------------------------
 !Decide what kind of simulation this is:
if (loopinit .eq. 0) then
   !This is a new simulation starting from data in qq_init.r8

   !Start from t = 0:
  tinit=zero
  t=zero

   !Read in initial PV (into qs) and compute its domain average:
  open(11,file='qq_init.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,qs
  call average(qs,qavg)
   !qavg is preserved throughout the simulation
  close(11)

   !Compute contour interval for PV:
  qqmax=maxval(qs)
  qqmin=minval(qs)
  dq=(qqmax-qqmin)/dble(ncontq)
   !Write information to log file:
  write(*,'(a,1x,f13.5,1x,f11.5,1x,f9.5)') &
          ' q_min, q_max, dq = ',qqmin,qqmax,dq
   !dq is preserved throughout the simulation

   !Copy qs into qd for recontouring:
  qd=qs

   !Initially there are no contours; set relevant counters to zero:
  nq=0
  nptq=0

else
   !In this case, we are restarting from existing data.  First read in
   !the contours and residual PV (see common.f90):
  call readcont(loopinit,1)
   !This also defines the entire PV field qs for restart and defines
   !the initial time, t.  It sets the value of tinit = t for use below.

endif

!--------------------------------------------------------------
 !Open files (depending on type of simulation):
if ((loopinit .eq. 0) .or. replace) then
   !This is either a new simulation starting from t = 0 (loopinit = 0)
   !or one which starts from a specified time in a simulation already
   !performed, and we only want to use the data at this time for the
   !initial conditions (no previous data is kept from the old simulation):

  if (replace) then
     !Obtain the last time in ene.asc BEFORE tinit found from contour files:
    open(15,file='evolution/ene.asc',status='old')
    td=zero
    do
      iread=0
      read(15,*,iostat=iread) td
      if (iread .ne. 0) exit 
      if (td .ge. tinit) exit
    enddo

    backspace(15)
    backspace(15)
    read(15,*) td
    write(*,*) ' Read ene.asc up to time t = ',td
    close(15)
     !Initialise counter and offset for writing direct access grid files:
    igrec=0
    igoff=1
  else
    td=zero
    igrec=0
    igoff=1
  endif

   !Set the correct initial time based on last entry in ene.asc:
  tinit=td
  t=tinit

   !Open all plain text diagnostic files (ending in ".asc"):
  open(12,file='evolution/monitor.asc',status='replace')
  open(14,file='contours/complexity.asc',status='replace')
  open(15,file='evolution/ene.asc',status='replace')
   !Open files for coarse grid saves (unformatted files ending in ".r4"):
  open(31,file='evolution/qq.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(32,file='evolution/zz.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
   !Open file for 1d PV spectra (ascii file):
  open(51,file='spectra/qspec.asc',status='replace')
   !Open file for contour writes (synopsis file):
  open(80,file='contours/qqsynopsis.asc',status='replace')

   !Initialise counter for writing direct access contour files:
  icrec=1

else
   !Here, we continue the simulation from time t = tinit, reading 
   !from the data in the cont subdirectory.  All previous data is
   !kept; we append more data on to existing files after making 
   !sure that they are synchronised.  Note: the base time for the
   !simulation is read from the first line of ene.asc.
  open(12,file='evolution/monitor.asc',status='old')
  td=zero
  do
    iread=0
    read(12,*,iostat=iread) td
    if (iread .ne. 0) exit 
    if (td .ge. tinit) exit
  enddo
  backspace(12)
  if (iread .eq. 0) then
    write(*,*) ' Read monitor.asc up to time before t = ',td
  else
    write(*,*) ' Read monitor.asc up to time t = ',td
  endif

  open(14,file='contours/complexity.asc',status='old')
  td=zero
  do
    iread=0
    read(14,*,iostat=iread) td
    if (iread .ne. 0) exit 
    if (td .ge. tinit) exit
  enddo
  backspace(14)
  if (iread .eq. 0) then
    write(*,*) ' Read complexity.asc up to time before t = ',td
  else
    write(*,*) ' Read complexity.asc up to time t = ',td
  endif

  open(15,file='evolution/ene.asc',status='old')
   !Obtain initial "base" time for simulation:
  read(15,*) tbase
  backspace(15)

  td=zero
  nreads=0
  do
    iread=0
    read(15,*,iostat=iread) td
    if (iread .ne. 0) exit 
    if (td .ge. tinit) exit
    nreads=nreads+1
  enddo
  backspace(15)
  if (iread .eq. 0) then
    write(*,*) ' Read ene.asc up to time before t = ',td
  else
    write(*,*) ' Read ene.asc up to time t = ',td
  endif

   !Open files for coarse grid saves:
  open(31,file='evolution/qq.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)
  open(32,file='evolution/zz.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)

   !Open file for 1d PV spectra:
  open(51,file='spectra/qspec.asc',status='old')
  do m=1,nreads
    read(51,*) td
    do k=1,kmaxred
      read(51,*)
    enddo
  enddo
  write(*,*) ' Read qspec.asc up to t = ',td

   !Open file for contour writes:
  open(80,file='contours/qqsynopsis.asc',status='old')
  do m=1,loopinit
    read(80,*)
  enddo

   !Initialise counters for writing direct access files:
  igrec=nreads
  igoff=1
  icrec=loopinit+1

  write(*,*) ' igrec = ',igrec
  write(*,*) ' igoff = ',igoff
  write(*,*) ' icrec = ',icrec

   !Set initial time to the base time in ene.asc:
  tinit=tbase
  write(*,*) ' tinit = ',tinit

endif

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect PV until next recontouring or end:
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
call recontour(qd)
write(*,'(a,i8,a,i9,a,f13.5)') ' nq = ',nq,'  nptq = ',nptq,'  t = ',t

return 
end subroutine

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(12) 
close(14)
close(15)
close(31)
close(32)
close(51)
close(80)

return
end subroutine

 !End main program
end program
!=======================================================================
