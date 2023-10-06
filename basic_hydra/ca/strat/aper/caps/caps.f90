!#####################################################################
!            The Combined Lagrangian Advection Method for
!        2D non-rotating Boussinesq flow in channel geometry
!#####################################################################

!        Code written by Stuart King & David Dritschel @ St Andrews
!                ***Version 1.0 completed 7 September 2012***

!        Principally adapted from bsl.f90; a fully SL solver for the problem 
!        described below.

!        This code solves: 
!            D zeta/ Dt = cos(th)*db/dx-sin(th)*db/dy = S_zz  
!               Db / Dt = 0
!               div (u) = 0
!        in the domain xmin < x < xmax ; ymin < y < ymax
!        (free slip boundary conditions, where th(t) is a simple
!         harmonic function oscillating between specified angles 
!         at a specified frequency).

!     This is treated using a semi-Lagrangian scheme for vorticity. 
!     Buoyancy is dealt with via a contour scheme. This allows us to  
!     compute b at each new timestep. This information can then be used
!     to accurately obtain the source term for the zeta equation. The
!     zeta equation can then be integrated forward a timestep using the 
!     computed source S_zz. 
!     Incompressibility is handled via the introduction of a stream-
!     function such that:
!             Lap(psi) = zeta
!     Determining psi (hence u,v) at each time step is termed the 
!     'inversion' problem.  Here this is done using a multi-grid
!     method, with psi = 0 boundary conditions.

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
do while (t .le. tfin)

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

implicit double precision(a-h,o-z)
implicit integer(i-n)

!-----------------------------------------------------------------
 !Read in initial vorticity (zs) and buoyancy (bb):
open(11,file='zz_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) dum,zs
close(11)
open(12,file='bb_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(12,rec=1) t,bb
close(12)

 !Compute contour interval for buoyancy:
bbmax=zero
bbmin=zero
do ix=0,nx
  do iy=0,ny
    bbmax=max(bbmax,bb(iy,ix))
    bbmin=min(bbmin,bb(iy,ix))
  enddo
enddo
bjump=(bbmax-bbmin)/dble(ncontb)

 !Write information to log file:
write(*,*)
write(*,'(a,3(1x,f9.5))') ' b_min, b_max, bjump = ',bbmin,bbmax,bjump

 !Copy zs into zd for recontouring:
do ix=0,nx
  do iy=0,ny
    zd(iy,ix)=zs(iy,ix)
  enddo
enddo

 !Compute contour interval for vorticity:
zzl1=zero
zzl2=zero
do ix=1,nxm1
  zzl1=zzl1+abs(zs(0,ix))+abs(zs(ny,ix))
  zzl2=zzl2+zs(0,ix)**2+zs(ny,ix)**2
enddo
do iy=1,nym1
  zzl1=zzl1+abs(zs(iy,0))+abs(zs(iy,nx))
  zzl2=zzl2+zs(iy,0)**2+zs(iy,nx)**2
enddo
zzl1=f12*zzl1+f14*(abs(zs(0,0))+abs(zs(ny,0))+abs(zs(0,nx))+abs(zs(ny,nx)))
zzl2=f12*zzl2+f14*(zs(0,0)**2+zs(ny,0)**2+zs(0,nx)**2+zs(ny,nx)**2)

do ix=1,nxm1
  do iy=1,nym1
    zzl1=zzl1+abs(zs(iy,ix))
    zzl2=zzl2+zs(iy,ix)**2
  enddo
enddo

zzl1=garea*zzl1
zzl2=garea*zzl2
 !Note: garea is the grid box area, glx*gly

if (zzl1 .gt. zero) then 
  zjump=(zzl2/zzl1)/dble(ncontz)
else
  zjump=zero
endif

!------------------------------------------------------------
 !Initially there are no contours:
nb=0
nptb=0
nz=0
nptz=0

 !Initialise time step so that subroutine adapt chooses a suitable one:
dt=zero

 !Set final time for simulation end:
itime=int((t+small)/tgsave)
tgrid=tgsave*dble(itime)
tfin=tgrid+tsim

 !Set flag for initialising reference state potential energy:
if (osci) then
   !Domain is oscillating; we don't compute a reference PE
  iene=1
  arg=omegag*t
  theta=theini+atilt*cos(arg)+btilt*sin(arg)
  ctheta=cos(theta)
  stheta=sin(theta)  
else
   !Domain is static 
  iene=0
  theta=theini
  ctheta=cos(theta)
  stheta=sin(theta)
endif

!--------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral(bbmax-bbmin)
 !Initialise constants and arrays for contour advection:
call init_contours

!--------------------------------------
 !Open all diagnostic plain text  files:     
open(12,file='monitor.asc',status='unknown')
open(13,file='norms.asc',status='unknown')
open(14,file='complexity.asc',status='unknown')
open(15,file='ene.asc',status='unknown')
 !Open file for 1d vorticity spectrum:
open(50,file='zspec.asc',status='unknown')
 !Open files for vorticity and buoyancy coarse
 !grid direct writes:
open(31,file='zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(32,file='bb.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
 !Open files for contour writes:
open(80,file='cont/zzsynopsis.asc',status='unknown')
open(83,file='cont/zzresi.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(90,file='cont/bbsynopsis.asc',status='unknown')

 !Initialise counter for writing direct files to the correct counter:
igrids=0

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect buoyancy & vorticity until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine

!=======================================================================

subroutine recont

use congen

implicit none

!Obtain new buoyancy contours:
write(*,*) 'Recontouring buoyancy ...'
call recontour(bb,xb,yb,bjump,bavg,nextb,indb,npb,i1b,i2b,nb,nptb,0)
write(*,'(a,i8,a,i9,a,f9.5)') '   nb = ',nb,'   nptb = ',nptb,'   db = ',bjump

!Obtain new vorticity contours:
write(*,*) 'Recontouring vorticity ...'
call recontour(zd,xz,yz,zjump,zavg,nextz,indz,npz,i1z,i2z,nz,nptz,1)
write(*,'(a,i8,a,i9,a,f9.5)') '   nz = ',nz,'   nptz = ',nptz,'   dz = ',zjump

return 
end subroutine

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(12) 
close(13)
close(14)
close(15)
close(31)
close(32)
close(50)
close(80)
close(83)
close(90)

return
end subroutine

 !End main program
end program
!=======================================================================

