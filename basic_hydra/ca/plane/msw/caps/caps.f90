!#########################################################################
!         The Doubly-Periodic Single-Layer Magneto-Hydrostatic
!          Shallow-Water Combined Lagrangian Advection Method
!#########################################################################

!       Code developed in July 2021 by D G Dritschel @ St Andrews.

!       This code simulates the unforced Magneto-Hydrostatic Shallow-Water
!       Equations (MSWE) in variables (q_l,delta,gamma_l,b_x,b_y,b_z),
!       where q_l is the linearised potential vorticity (zeta - f*h_tilde),
!       delta is the velocity divergence, gamma_l is the linearised
!       acceleration divergence (f*zeta - c^2*Lap{h_tilde}), b_x & b_y
!       are the horizontal components of the magnetic field, and b_z is
!       the vertical component at the bottom z = 0.  Here h_tilde is the
!       dimensionless height anomaly, c^2 = g*H (H is the mean depth)
!       and zeta is the relative vorticity (dv/dx - du/dy).

!       Contour advection and generation are done internally now.  For
!       details of the method, see Dritschel & Fontane, J. Comput. Phys.
!       229, pp. 5408--5417 (2010).  Only q_l is represented by contours.

!       The full algorithm consists of the following modules:

!       caps.f90      : This source - main program loop, repeats successive 
!                       calls to evolve fields and recontour;
!       parameters.f90: User defined parameters for a simulation;
!       constants.f90 : Fixed constants used throughout the other modules;
!       common.f90    : Common data preserved throughout simulation 
!                       (through recontouring--evolution cycle);
!       spectral.f90  : Fourier transform common storage and routines;
!       contours.f90  : Contour advection common storage and routines;
!       congen.f90    : Source code for contour-to-grid conversion;
!       evolution.f90 : Main time evolution module - advects gridded 
!                       fields using a PS method along with contours.
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
double precision:: qmin,ql1,ql2
integer:: ix,iy

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!----------------------------------------------------------------------
 !Read in q_l as qr and convert to spectral space as qs:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,qr
close(11)
 !Note: qr must have a zero domain average, mathematically.

 !Choose PV contour interval based on <q^2>/<|q|> for |q| > q_rms:
qmin=sqrt(dsumi*sum(qr**2))
ql1=zero
ql2=zero
do ix=1,ng
  do iy=1,ng
    if (abs(qr(iy,ix)) .gt. qmin) then
      ql1=ql1+abs(qr(iy,ix))
      ql2=ql2+qr(iy,ix)**2
    endif
  enddo
enddo
qjump=ql2/(ql1*dble(ncont))

 !Convert qr to spectral space as qs (zz is overwritten):
zz=qr
call ptospc(ng,ng,zz,qs,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
qs(1,1)=zero

!----------------------------------------------------------------------
 !Read in delta and convert to spectral space as ds:
open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,ds,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
ds(1,1)=zero

!----------------------------------------------------------------------
 !Read in gamma_l and convert to spectral space as gs:
open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,gs,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
gs(1,1)=zero

!----------------------------------------------------------------------
 !Read in x component of B and convert to spectral space as bxs:
open(11,file='bx_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,bxs,xfactors,yfactors,xtrig,ytrig)

!----------------------------------------------------------------------
 !Read in y component of B and convert to spectral space as bys:
open(11,file='by_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,bys,xfactors,yfactors,xtrig,ytrig)

!----------------------------------------------------------------------
 !Read in z component of B(x,y,0) and convert to spectral space as bzs:
open(11,file='bz_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,bzs,xfactors,yfactors,xtrig,ytrig)

!----------------------------------------------------------------------
 !Spectrally-truncate all fields for use in de-aliasing:
qs=filt*qs
ds=filt*ds
gs=filt*gs
bxs=filt*bxs
bys=filt*bys
bzs=filt*bzs

 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv)
 !and relative vorticity (zz):
call main_invert(qs,ds,gs,hh,uu,vv,zz)
 !Note: qs, ds & gs are in spectral space while 
 !      hh, uu, vv and zz are in physical space.

!----------------------------------------------------------------------
 !Initially there are no contours (they are built by congen from qr):
nq=0
nptq=0

 !Initialise counters for saving gridded and contour data:
igrids=0
iconts=0

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='contours/complexity.asc',status='replace')
open(15,file='evolution/ecomp.asc',status='replace')
open(16,file='evolution/ro-fr-hm.asc',status='replace')
open(17,file='evolution/norms.asc',status='replace')
open(18,file='evolution/magdiss.asc',status='replace')

 !Open files for coarse grid saves:
open(32,file='evolution/dd.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(33,file='evolution/gg.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(34,file='evolution/hh.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(35,file='evolution/zz.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(36,file='evolution/jz.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(37,file='evolution/bp.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(38,file='evolution/bz.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)

 !Open files for 1d spectra:
open(51,file='spectra/zspec.asc',status='replace') ! zeta
open(52,file='spectra/dspec.asc',status='replace') ! delta
open(53,file='spectra/gspec.asc',status='replace') ! gamma_l
open(54,file='spectra/hspec.asc',status='replace') ! h_tilde
open(55,file='spectra/jspec.asc',status='replace') ! j_z
open(56,file='spectra/cspec.asc',status='replace') ! b_z

 !Open files for PV contour writes:
open(80,file='contours/qqsynopsis.asc',status='replace')
open(83,file='contours/qqresi.r8',form='unformatted',access='direct', &
                                status='replace',recl=nbytes)

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
call recontour(qr)
write(*,'(a,i8,a,i9)') '   nq = ',nq,'   nptq = ',nptq

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
close(18)

close(32)
close(33)
close(34)
close(35)
close(36)
close(37)
close(38)

close(51)
close(52)
close(53)
close(54)
close(55)
close(56)

close(80)
close(83)

return
end subroutine finalise

 !End main program
end program caps
!=======================================================================
