!#########################################################################
!        The Doubly-Periodic Single-Layer Quasi-Geostrophic (QG)
!             Incompressible Magneto-Hydrodynamical (IMHD)
!              Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!        Code written by David Dritschel on 6 Dec 2013 @ St Andrews
!        Revised to restructure data files and to allow for spectral
!        forcing by DGD in late February 2021.

!        This code simulates the following system of equations:

!              Dq/Dt = S_u + S_b + F_q
!              DA/Dt = -B_0*dpsi/dx + eta*Lap{A} + F_A

!        where q is the QG potential vorticity (PV) and A is the magnetic
!        potential, related to the magnetic field B via

!              curl(A*e_z) = (A_y, -A_x, 0).

!        Above, B_0 is a mean magnetic field in the x direction, eta is
!        the magnetic diffusivity, while F_q & F_A allow for stochastic
!        forcing of vorticity and current density j = -Lap(A).

!        ==================================================
!        ====>  The magnetic Prandtl number is zero.  <====
!        ==================================================

!        The PV source terms are

!             S_u = -r*zeta + psi/(tau*L_D^2)       and
!             S_b = J(A,Lap{A}) - B_0*Lap{dA/dx}

!        The velocity field (u,v) is found by inverting

!             Lap{psi} - psi/L_D^2 = q - beta*y
!             u = -dpsi/dy ; v = dpsi/dx

!        where
!             r      is the Ekman damping rate
!             tau    is the thermal damping rate
!             L_D    is the Rossby deformation length
!             beta   is the planetary vorticity gradient

!        The forcing (if present) is Markovian, time correlated, and
!        maintained at a prescribed rms input rate.  The spectrum of
!        the forcing is c*k^{2p-1}*exp(-2*(k/k_0)^2), where p and
!        k_0 are specified in parameters.f90, and c is determined by
!        the specified enstrophy (variance) input rate (esrz or esrj).

!        We split the PV evolution equation into *three* equations,

!             Dq_c/Dt = 0  :  contour advection
!             Dq_s/Dt = 0  :  pseudo-spectral
!             Dq_d/Dt = S  :  pseudo-spectral

!        such that q is a weighted sum of these,

!             q = F(q_s) + (1-F)*(q_c) + q_d

!        where F is a low-pass filter defined in spectral.f90 and 1-F is
!        a complementary high pass filter (see Dritschel & Fontane, JCP,
!        2010).

!        Hence advection at large to intermediate scales is controlled 
!        by the pseudo-spectral method, whereas advection at intermediate
!        to small scales is controlled by contour advection (where it 
!        is most accurate).  The source terms S_u & S_b are handled 
!        entirely by the pseudo-spectral method (through q_d).

!        At the beginning of each time step, q_s is replaced by q,
!        while q_d is replaced by (1-F)(q-q_c).  q_c remains as
!        contours for a period of time determined by twistmax below.
!        After this period, q is obtained on an ultra-fine grid and
!        re-contoured so that the accumulated forcing in q_d is
!        given to the contours in q_c (to the extent possible).

!        Time integration uses the 4th-order Runge-Kutta method, with
!        diffusion built in using spectral integrating factors.

!        The full algorithm consists of the following modules:

!        casl.f90      : This source - main program loop, repeats
!                        successive calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        common.f90    : Common data preserved throughout simulation 
!                        (through the recontouring-evolution cycle);
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded 
!                        fields using a PS method along with contours.
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

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

 !Local variables:
double precision:: wka(ny,nx),qmin,ql1,ql2
integer:: ix,iy,itime

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!----------------------------------------------------------------------
 !Read in PV, possibly subtract beta*y, and convert to spectral space:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qr
close(11)

if (beffect) then
   !Subtract beta*y (in bety below) to define the PV anomaly:
  do ix=1,nx
    qr(:,ix)=qr(:,ix)-bety
  enddo
   !Choose an integral number of PV jumps (required in recontour):
  qjump=beta*hly/dble(ncontq)
else
   !Choose contour interval based on <q^2>/<|q|> for |q| > q_rms:
  qmin=sqrt(dsumi*sum(qr**2))
  ql1=zero
  ql2=zero
  do ix=1,nx
    do iy=1,ny
      if (abs(qr(iy,ix)) .gt. qmin) then
        ql1=ql1+abs(qr(iy,ix))
        ql2=ql2+qr(iy,ix)**2
      endif
    enddo
  enddo
  qjump=ql2/(ql1*dble(ncontq))
endif

 !Copy qr into wka before taking Fourier transform (wka is overwritten):
wka=qr !qr is needed below for recontouring
call ptospc(nx,ny,wka,qs,xfactors,yfactors,xtrig,ytrig)
 !Ensure mean value is zero:
qs(1,1)=zero

!----------------------------------------------------------------------
 !Read in magnetic potential, A, and convert to spectral space:
open(11,file='aa_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,wka
close(11)
call ptospc(nx,ny,wka,aa,xfactors,yfactors,xtrig,ytrig)
 !Ensure mean value is zero:
aa(1,1)=zero

!----------------------------------------------------------------------
 !Initially there are no contours (they are built by congen from qr):
nq=0
nptq=0

 !Initialise counter for saving gridded data:
igrids=0

 !Initialise vorticity and/or magnetic potential forcing (if present):
if (zforcing) call ranspec(dzdt,esrz,powz,k0z)
if (aforcing) call ranspec(dadt,esra,powa,k0a)

!----------------------------------------------------------------------
 !Open all plain text diagnostic files:
open(14,file='contours/complexity.asc',status='replace')
open(15,file='evolution/ecomp.asc',status='replace')
open(16,file='evolution/norms.asc',status='replace')
open(17,file='evolution/monitor.asc',status='replace')
if (tracer) open(24,file='evolution/circulation.asc',status='replace')

 !Open files for coarse grid saves:
open(31,file='evolution/qq.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(32,file='evolution/jj.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(33,file='evolution/aa.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)

 !Open files for 1d PV & current density spectra:
open(51,file='spectra/qspec.asc',status='replace')
open(52,file='spectra/jspec.asc',status='replace')

 !Open files for contour writes:
open(80,file='contours/qqsynopsis.asc',status='replace')
open(83,file='contours/qqresi.r4',form='unformatted',access='direct', &
                                status='replace',recl=nbytes)

return
end subroutine initialise

!=======================================================================

subroutine evolve

use evolution

implicit none

 !Advect q & A until next recontouring or end:
write(*,*) ' Evolving contours and fields ...'
call advect

return 
end subroutine evolve

!=======================================================================

subroutine recont

use congen

implicit none

!Obtain new PV contours:
write(*,*) ' Recontouring PV ...'
call recontour
write(*,'(a,i8,a,i9,a,f9.5)') '   nq = ',nq,'   nptq = ',nptq,'   dq = ',qjump

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
close(17)
if (tracer) close(24)
close(31)
close(32)
close(33)
close(51)
close(52)
close(80)
close(83)

return
end subroutine finalise

 !End main program
end program caps
!=======================================================================
