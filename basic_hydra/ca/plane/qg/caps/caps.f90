!#########################################################################
!                    The Doubly-Periodic Single-Layer 
!        Quasi-Geostrophic Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!        Code adapted from the code in imhd on 17 Jan 2014 @ St Andrews

!          This code simulates the following system of equations:

!             Dq/Dt = S 

!          where q is the QG PV. The PV source term is

!             S = -r*zeta + (psi-psi_eq)/(tau*L_D^2) + F

!          The velocity field (u,v) is found by inverting 

!             Lap{psi} - psi/L_D^2 = q - beta*y   
!             u = -dpsi/dy ; v = dpsi/dx                   

!          where
!             r      is the Ekman damping rate
!             tau    is the thermal damping rate
!             L_D    is the Rossby deformation length
!             psi_eq is the thermal equilibrium streamfunction
!             F      is stochastic forcing (by vortex injection)
!             beta   is the planetary vorticity gradient

!          We split the PV evolution equation into *three* equations,

!             Dq_c/Dt = 0  :  contour advection
!             Dq_s/Dt = 0  :  pseudo-spectral
!             Dq_d/Dt = S  :  pseudo-spectral

!          such that q is a weighted sum of these,

!             q = F(q_s) + (1-F)*(q_c) + q_d

!          where F is a low-pass filter defined in spectral.f90 and 1-F is
!          a complementary high pass filter (see Dritschel & Fontane, JCP,
!          2010).

!          Hence advection at large to intermediate scales is controlled 
!          by the pseudo-spectral method, whereas advection at intermediate
!          to small scales is controlled by contour advection (where it 
!          is most accurate).  The source term S is handled entirely by
!          the pseudo-spectral method.

!          At the beginning of each time step, q_s is replaced by q,
!          while q_d is replaced by (1-F)(q-q_c).  q_c remains as
!          contours for a period of time determined by twistmax below.
!          After this period, q is obtained on an ultra-fine grid and
!          re-contoured so that the accumulated forcing in q_d is
!          given to the contours in q_c (to the extent possible).

!     The full algorithm consists of the following modules:
!        casl.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        generic.f90   : Generic service routines for CASL;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded 
!                        fields using a PS method along with contours.
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

   !Advect PV and other fields until next recontouring or end:
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

implicit none

 !Local variables:
double precision:: ff(ny,nx)
double precision:: uni,tbase,td
integer:: i,k,m,ix,iy,iread,nreads

!------------------------------------------------------------------
 !Allocate memory for thermal equilibrium streamfunction if heating
 !is present:
if (heating) allocate(ppeq(nx,ny))

!--------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

 !Read in thermal equilibrium streamfunction (if present) and
 !convert to spectral coefficients (ppeq):
if (heating) then
  open(12,file='psieq.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) td,ff
  close(12)
  call ptospc(nx,ny,ff,ppeq,xfactors,yfactors,xtrig,ytrig)
endif

!-----------------------------------------------------------------
 !Decide what kind of simulation this is:
if (loopinit .eq. 0) then
   !This is a new simulation starting from data in qq_init.r8

   !Start from t = 0:
  tinit=zero

   !Read in full PV, possibly subtract beta*y, and FFT to initialise qs:
  open(11,file='qq_init.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,qr
  close(11)

  if (beffect) then
     !Subtract beta*y to define the PV anomaly:
    do ix=1,nx
      do iy=1,ny
        qr(iy,ix)=qr(iy,ix)-bety(iy)
      enddo
    enddo
     !Choose an integral number of PV jumps (required in recontour):
    qjump=beta*elly/dble(ncontq)
    write(*,*)
    write(*,'(a,f14.8,a,f12.8)') ' beta = ',beta,'   qjump = ',qjump

  else
     !Choose contour interval based on range of PV values:
    call contint(qr,ncontq,qjump)
    write(*,*)
    write(*,'(a,1x,f13.8)') ' qjump = ',qjump
  endif

   !Copy PV anomaly qr into ff before taking Fourier transform 
   !(ff is overwritten):
  ff=qr
   !Convert PV anomaly to spectral space as qs:
  call ptospc(nx,ny,ff,qs,xfactors,yfactors,xtrig,ytrig)

   !Initially there are no contours:
  nq=0
  nptq=0

else
   !In this case, we are restarting from existing data.  First read in
   !the contours and residual PV (see common.f90 for readcont):
  call readcont(loopinit)
   !This also defines the entire PV field qs for restart and defines
   !the initial time, t.  It sets the value of tinit = t for use below.
endif

!--------------------------------------------------------------------
 !If stochastic forcing is used, initialise vortex forcing variables:
if (stoch) then
  if (ivor .eq. 1) then
     !Point vortices of mean (grid) vorticity of +/-vorvor are
     !added at an average rate of dnvor vortices per unit time:
    dnvor=two*esr*(three*pi/vorvor)**2/garea
  else
     !Point vortex dipoles (concentrated to points) are added
     !at an average rate of dnvor dipoles per unit time, with 
     !a maximum absolute grid vorticity of vorvor:
    dnvor=six*esr*(pi/vorvor)**2/garea
  endif
   !Above, esr is the enstrophy input rate.

   !Initialize random # generator on first call:
  do i=1,iseed
    uni=rand(0)
  enddo

   !Initialise total added vortices:
  totnvor=tinit*dnvor
endif

 !Initialise time step so that subroutine adapt chooses a suitable one:
dt=zero

!--------------------------------------------------------------
 !Open files (depending on type of simulation):
if (loopinit .eq. 0) then
   !This is a new simulation starting from t = 0
  tinit=zero
  t=zero

   !Initialise counter for writing direct access grid files:
  igrec=0

   !Initialise counter for writing direct access contour files:
  icrec=1

   !Open all plain text diagnostic files:
  open(14,file='complexity.asc',status='replace')
  open(15,file='ene.asc',status='replace')
  open(16,file='norms.asc',status='replace')
  open(17,file='monitor.asc',status='replace')

   !Open file for 1d PV spectra:
  open(51,file='spectra.asc',status='replace')

   !Open file for coarse grid saves:
  open(31,file='qq.r4',form='unformatted',access='direct', &
                   & status='replace',recl=nbytes)

   !Open files for contour writes:
  open(80,file='cont/qqsynopsis.asc',status='replace')
  open(83,file='cont/qqresi.r4',form='unformatted',access='direct', &
                            & status='replace',recl=nbytes)

else
   !Here, we continue the simulation from time t = tinit, reading 
   !from the data in the cont subdirectory.  All previous data is
   !kept; we append more data on to existing files after making 
   !sure that they are synchronised.  Note: the base time for the
   !simulation is read from the first line of ene.asc.
  open(17,file='monitor.asc',status='old')
  td=zero
  do
    iread=0
    read(17,*,iostat=iread) td
    if (iread .ne. 0) exit 
    if (td .ge. tinit) exit
  enddo
  backspace(17)
  if (iread .eq. 0) then
    write(*,*) ' Read monitor.asc up to time before t = ',td
  else
    write(*,*) ' Read monitor.asc up to time t = ',td
  endif

  open(16,file='norms.asc',status='old')
  td=zero
  do
    iread=0
    read(16,*,iostat=iread) td
    if (iread .ne. 0) exit 
    if (td .ge. tinit) exit
  enddo
  backspace(16)
  if (iread .eq. 0) then
    write(*,*) ' Read norms.asc up to time before t = ',td
  else
    write(*,*) ' Read norms.asc up to time t = ',td
  endif

  open(14,file='complexity.asc',status='old')
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

  open(15,file='ene.asc',status='old')
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
   !Initialise contour for writing data every tgsave time units:
  if (iread .eq. 0) then
    write(*,*) ' Read ene.asc up to time before t = ',td
    igrec=nreads+1
  else
    write(*,*) ' Read ene.asc up to time t = ',td
    igrec=nreads
  endif

   !Open file for coarse grid saves:
  open(31,file='qq.r4',form='unformatted',access='direct', &
                   & status='old',recl=nbytes)

   !Open file for 1d PV spectra:
  open(51,file='spectra.asc',status='old')
  do m=1,nreads
    read(51,*) td
    do k=1,kmaxred
      read(51,*)
    enddo
  enddo
  write(*,*) ' Read spectra.asc up to t = ',td

   !Open file for contour writes:
  open(80,file='cont/qqsynopsis.asc',status='old')
  do m=1,loopinit
    read(80,*)
  enddo
   !Open file containing residual PV:
  open(83,file='cont/qqresi.r4',form='unformatted',access='direct', &
                            & status='old',recl=nbytes)

   !Initialise counters for writing direct access contour files:
  icrec=loopinit+1

  write(*,*) ' igrec = ',igrec
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

if (forcing) call contint(qr,ncontq,qjump)

call recontour(qr)
write(*,'(a,i9,a,i10,a,f9.5)') '   nq = ',nq,'   nptq = ',nptq,'   dq = ',qjump

return 
end subroutine

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(14)
close(15)
close(16)
close(17)
close(31)
close(51)
close(80)
close(83)

return
end subroutine

 !End main program
end program
!=======================================================================
