!#########################################################################
!            The Combined Lagrangian Advection Method for
!      2D non-rotating Boussinesq flow in a periodic channel geometry
!#########################################################################

!        Code written by Stuart King & David Dritschel @ St Andrews
!                ***Version 1.0 completed 18 December 2012***
!                ***Updated by dgd on 15 January 2022***

!        This code evolves vorticity zeta and buoyancy b using
!               Dzeta/Dt = db/dx = S_zz  
!               Db/Dt = 0
!               div(u,v) = 0
!        in the domain xmin < x < xmax ; ymin < y < ymax
!        (free slip boundary conditions in y, periodic in x).

!     Vorticity is handled by a combination of a pseudo-spectral scheme
!     and contour advection (CLAM - see Dritschel & Fontane, JCP 2010).
!     Buoyancy is handled purely by contour advection as it is exactly
!     conserved here.  The vorticity source S_zz is obtained directly
!     from the buoyancy contours, effectively integrating the delta-
!     function contributions associated with each buoyancy contour
!     Details can be found in section 3.3 of Carr et al, JFM 2017.

!     Incompressibility is exactly satisfied using a streamfunction
!     psi satisfying
!             Lap(psi) = zeta
!     Determining psi (hence u,v) at each time step is termed the 
!     'inversion' problem.  Here this is done using a spectral method
!     including exact solutions of the homogeneous problem Lap(psi) = 0
!     to satisfy the boundary conditions (v = 0).

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
!        evolution.f90 : Main time evolution module - advects gridded fields 
!                        using SL method along with contours.
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

   !Obtain new buoyancy & vorticity contours:
  call recont

   !Advect buoyancy & vorticity until next recontouring or end:
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

double precision:: bbmax,bbmin,zzmin,zzl1,zzl2
integer:: ix,iy

!-----------------------------------------------------------------
 !Read in initial vorticity (zz) and buoyancy (bb):
open(11,file='zz_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)

open(11,file='bb_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,bb
close(11)

!-----------------------------------------------------------------
 !Compute contour interval for buoyancy:
bbmax=maxval(bb)
bbmin=minval(bb)
bjump=(bbmax-bbmin)/dble(ncontb)

 !Write information to log file:
write(*,*)
write(*,'(a,3(1x,1p,e14.7))') ' b_min, b_max, bjump = ',bbmin,bbmax,bjump

!-----------------------------------------------------------------
 !Compute domain-average buoyancy and vorticity (invariant):
bavg=dsumi*(f12*sum(bb(0,:)+bb(ny,:))+sum(bb(1:nym1,:)))
zavg=dsumi*(f12*sum(zz(0,:)+zz(ny,:))+sum(zz(1:nym1,:)))
 !Note: dsumi = 1/(nx*ny)

!-----------------------------------------------------------------
 !Choose vorticity contour interval based on <zeta^2>/<|zeta|>
 !for |zeta| > zeta_rms:
zzmin=sqrt(dsumi*(f12*sum(zz(0,:)**2+zz(ny,:)**2)+sum(zz(1:nym1,:)**2)))
 !dsumi = 1/(nx*ny)
if (zzmin .gt. zero) then 
  zzl1=zero
  zzl2=zero
  do ix=0,nxm1
    if (abs(zz(0,ix)) .gt. zzmin) then
      zzl1=zzl1+abs(zz(0,ix))
      zzl2=zzl2+zz(0,ix)**2
    endif
    if (abs(zz(ny,ix)) .gt. zzmin) then
      zzl1=zzl1+abs(zz(ny,ix))
      zzl2=zzl2+zz(ny,ix)**2
    endif
  enddo
  zzl1=f12*zzl1
  zzl2=f12*zzl2
  do ix=0,nxm1
    do iy=1,nym1
      if (abs(zz(iy,ix)) .gt. zzmin) then
        zzl1=zzl1+abs(zz(iy,ix))
        zzl2=zzl2+zz(iy,ix)**2
      endif
    enddo
  enddo
  zjump=zzl2/(zzl1*dble(ncontz))
else
  zjump=zero
endif

!-----------------------------------------------------------------
 !Initially there are no contours:
nb=0
nptb=0
nz=0
nptz=0

 !Set flag to compute reference potential energy (peref) in contours.f90:
iene=0

!-----------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral

 !Initialise constants and arrays for contour advection:
call init_contours

 !Convert vorticity to spectral space as zs:
uu=zz !Here uu is used temporarily
call ptospc_fc(nx,ny,uu,zs,xfactors,yfactors,xtrig,ytrig)
zs(0,0)=zero !Set mean value to zero (the mean is not needed for zs)

!-----------------------------------------------------------------
 !Open all plain text diagnostic files:
open(14,file='evolution/complexity.asc',status='unknown')
open(21,file='evolution/ecomp.asc',status='replace')
open(22,file='evolution/monitor.asc',status='replace')
open(23,file='evolution/vorticity.asc',status='replace')
 !Open files for coarse grid saves:
open(31,file='evolution/zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(32,file='evolution/bb.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
 !Open file for 1d vorticity & buoyancy spectra:
open(51,file='spectra/zspec.asc',status='unknown')
open(52,file='spectra/bspec.asc',status='unknown')
 !Open files for contour writes:
open(80,file='cont/zzsynopsis.asc',status='unknown')
open(83,file='cont/zzresi.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(90,file='cont/bbsynopsis.asc',status='unknown')

!-----------------------------------------------------------------
 !Initialise counters for saving gridded and contour data:
igrids=0
iconts=0

return
end subroutine initialise

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect buoyancy & vorticity until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine evolve

!=======================================================================

subroutine recont

use congen

implicit none

!Obtain new buoyancy contours:
write(*,*) 'Recontouring buoyancy ...'
call recontour(bb,xb,yb,bjump,bavg,nextb,indb,npb,i1b,i2b,nb,nptb,0)
write(*,'(a,i8,a,i9,a,1p,e14.7)') '   nb = ',nb,'   nptb = ',nptb,'   db = ',bjump

!Obtain new vorticity contours:
write(*,*) 'Recontouring vorticity ...'
call recontour(zz,xz,yz,zjump,zavg,nextz,indz,npz,i1z,i2z,nz,nptz,1)
write(*,'(a,i8,a,i9,a,1p,e14.7)') '   nz = ',nz,'   nptz = ',nptz,'   dz = ',zjump

return
end subroutine recont

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(14)
close(21)
close(22)
close(23)
close(31)
close(32)
close(51)
close(52)
close(80)
close(83)
close(90)

return
end subroutine finalise

 !End main program
end program caps
!=======================================================================
