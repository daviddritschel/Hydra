!###########################################################################
!    The ellipsoidal barotropic Contour-Advective Semi-Lagrangian (CASL) 
!    algorithm for 2D flow simulations on an ellipsoid of revolution.

!    Adapted from the spherical version in January 2015 by D.G. Dritschel
!###########################################################################

!       Uses a 4th-order Runge-Kutta scheme with an adapted time step
!       together with 4th-order compact differencing in latitude.

!       Input data files:
!       =================
!        cont/synopsis.asc  initial time, number of contours and nodes
!         "   index000      contour counters, etc at this time
!         "   nodes000      contour nodes at this time

!       Output data files:
!       =================
!       Every nsteps time steps (or every "period"):
!        monitor.asc        time, energy, enstrophy & angular momentum
!        zz.r4              gridded relative vorticity field
!        cont/synopsis.asc  current time, number of contours and nodes
!         "   indexnnn      contour counters, etc at period "nnn" 
!         "   nodesnnn      contour nodes at period "nnn"

!       Every contour regularisation (surgery):
!        complexity.asc     number of contours, nodes and time

!==========================================================================

!     The full algorithm consists of the following modules:
!        casl.f90      : This source - main program loop, advects contours
!                        and regularises them by surgery; writes data;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!----------------------------------------------------------------------------
program casl

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use spectral
use contours

implicit none

 !Local variables:
integer,parameter:: nreg=8
 !nreg: every nreg time steps, surgery is performed (Fontane & D, JCP 2009)
integer:: ireg, loop, istep

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

 !Initialise (local) counter for surgery:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do loop=1,nperiod
  do istep=1,nsteps

     !Evolve contours from t to t + dt:
    call advance

    ireg=ireg+1
    if (ireg .gt. nreg) then
       !Regularise the contours (surgery and node redistribution):
      call surgery
       !Record contour complexity:
      write(30,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt
       !Re-set counter:
      ireg=0
    endif
  enddo  

   !Save data periodically:
  call writedata(loop)
enddo
 !End time loop 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Close all files:
call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

integer:: i,j

!--------------------------------------------------------------------
 !Initialise modules:
call init_spectral
call init_contours

 !Read initial contours:
open(80,file='cont/synopsis.asc',status='old')
read(80,*) t,n,npt

open(81,file='cont/index000',form='unformatted',status='old')
read(81) np(1:n),i1(1:n),ind(1:n)
close(81)

open(82,file='cont/nodes000',form='unformatted',status='old')
read(82) x(1:npt),y(1:npt),z(1:npt)
close(82)

 !Define contour end points (i2) together with the "next" array:
do j=1,n
  i2(j)=i1(j)+np(j)-1
  do i=i1(j),i2(j)-1
    next(i)=i+1
  enddo
  next(i2(j))=i1(j)
enddo 

!--------------------------------------------------------------------
 !Open various diagnostic files:

 !Energy & other diagnostics every time step:
open(20,file='monitor.asc',status='replace')

 !Contour complexity after surgery (n, npt & t):
open(30,file='complexity.asc',status='replace')

 !Relative vorticity every nsteps time steps:
open(40,file='zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !Initialise dump counter for writing to correct record:
idump=0

!--------------------------------------------------------------------
 !Obtain initial surface velocity field and write initial data:
call genvor
call invert
write(*,*) '       t            E              Z              M'
call writedata(0)

return
end subroutine

!=======================================================================

subroutine advance
! Integrates the equations of motion from time t to time t + dt.
! *** Uses the 4th-order Runge-Kutta method ***

implicit none

 !Local variables:
double precision:: u(npt),v(npt),w(npt)
double precision:: xi(npt),yi(npt),zi(npt)
double precision:: xf(npt),yf(npt),zf(npt)
double precision:: xnew,ynew,znew,fac
integer:: i

!---------------------------------------------------------------------
 !RK4 predictor step to time t0 + dt/2:

 !Obtain surface velocity field (uu,vv) from the contours:
call invert
 !Interpolate the velocity to the contour nodes (see contours.f90):
call velint(uu,vv,u,v,w)

do i=1,npt
  xi(i)=x(i)
  yi(i)=y(i)
  zi(i)=z(i)
  xnew=xi(i)+dt2*u(i)
  ynew=yi(i)+dt2*v(i)
  znew=zi(i)+dt2*w(i)
  fac=one/sqrt(xnew**2+ynew**2+znew**2)
  x(i)=fac*xnew
  y(i)=fac*ynew
  z(i)=fac*znew
  xf(i)=xi(i)+dt6*u(i)
  yf(i)=yi(i)+dt6*v(i)
  zf(i)=zi(i)+dt6*w(i)
enddo

!---------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt/2:
t=t+dt2
call invert
call velint(uu,vv,u,v,w)

do i=1,npt
  xnew=xi(i)+dt2*u(i)
  ynew=yi(i)+dt2*v(i)
  znew=zi(i)+dt2*w(i)
  fac=one/sqrt(xnew**2+ynew**2+znew**2)
  x(i)=fac*xnew
  y(i)=fac*ynew
  z(i)=fac*znew
  xf(i)=xf(i)+dt3*u(i)
  yf(i)=yf(i)+dt3*v(i)
  zf(i)=zf(i)+dt3*w(i)
enddo

!---------------------------------------------------------------------
 !RK4 predictor step at time t0 + dt:
call invert
call velint(uu,vv,u,v,w)

do i=1,npt
  xnew=xi(i)+dt*u(i)
  ynew=yi(i)+dt*v(i)
  znew=zi(i)+dt*w(i)
  fac=one/sqrt(xnew**2+ynew**2+znew**2)
  x(i)=fac*xnew
  y(i)=fac*ynew
  z(i)=fac*znew
  xf(i)=xf(i)+dt3*u(i)
  yf(i)=yf(i)+dt3*v(i)
  zf(i)=zf(i)+dt3*w(i)
enddo

!---------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt:
t=t+dt2
call invert
call velint(uu,vv,u,v,w)

do i=1,npt
  xnew=xf(i)+dt6*u(i)
  ynew=yf(i)+dt6*v(i)
  znew=zf(i)+dt6*w(i)
  fac=one/sqrt(xnew**2+ynew**2+znew**2)
  x(i)=fac*xnew
  y(i)=fac*ynew
  z(i)=fac*znew
enddo

return
end subroutine

!==========================================================================

subroutine invert
! Calculates the surface velocity field (uu,vv) given the vorticity contours.

implicit none

 !Local variables:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt)
integer:: i,j

!------------------------------------------------------------------------
 !Obtain the semi-spectral vorticity field:
call genvor

 !Invert Laplace's operator on the vorticity anomaly:
call laplinv(qq,wka,wkb)
 !Here the streamfunction psi is wka while wkb = tau*d(psi)/dlat.

 !Compute d(psi)/dlon = wkc:
call deriv(ng,nt,rk,wka,wkc)

 !Get physical space velocity:
call revfft(ng,nt,wkb,trig,factors)
call revfft(ng,nt,wkc,trig,factors)  

 !Copy into uu & vv which include polar points (see velint in contours.f90):
 !and correct meridional velocity by multiplying by tau/rho:
do i=1,nt
  do j=1,ng
    uu(j,i)=      -wkb(j,i)
    vv(j,i)=tdr(j)*wkc(j,i)
  enddo
enddo
 !Above, tdr = tau/rho.

return
end subroutine

!========================================================================

subroutine genvor
! Computes the vorticity field qq by performing a contour->grid conversion.
! ***qq is returned in semi-spectral space***

implicit none

 !Local variables:
integer:: i,j

!-----------------------------------------------------
 !Grid the vorticity contours as qq:
call con2grid(qq)

 !Subtract Coriolis frequency to define anomaly:
if (rotate) then
  do i=1,nt
    do j=1,ng
      qq(j,i)=qq(j,i)-cof(j)
    enddo
  enddo
endif

 !FFT qq in longitude:
call forfft(ng,nt,qq,trig,factors) 

return
end subroutine

!========================================================================

subroutine writedata(loop)
! This routine writes the current contours, the associated vorticity 
! field, and various diagnostics (energy, etc).

implicit none

 !Local variables:
double precision:: zz(ng,nt),wka(ng,nt),wkb(ng,nt),wkc(ng,nt)
double precision:: eke,ens,ang
real:: zzr4(ng,nt),tr4
integer:: i,j,m,loop
character(len=3):: ofile

!------------------------------------------------------------------
 !Get qq in physical space (as zz, the relative vorticity):
do m=1,nt
  do j=1,ng
    zz(j,m)=qq(j,m)
  enddo
enddo
call revfft(ng,nt,zz,trig,factors)

 !Compute energy, enstrophy and angular momentum per unit area:
do i=1,nt
  do j=1,ng
    wka(j,i)=rdt(j)*(uu(j,i)**2+tsqi(j)*vv(j,i)**2)
    wkb(j,i)=rdt(j)*zz(j,i)**2
    wkc(j,i)=rdt(j)*rho(j)*uu(j,i)
  enddo
enddo
 !Above, rdt = mu' = rho/tau while tsqi = 1/tau^2.

eke=zero
ens=zero
ang=zero
do i=1,nt
  eke=eke+f1112*(wka(1,i)+wka(ng,i))
  ens=ens+f1112*(wkb(1,i)+wkb(ng,i))
  ang=ang+f1112*(wkc(1,i)+wkc(ng,i))
  do j=2,ngm1
    eke=eke+wka(j,i)
    ens=ens+wkb(j,i)
    ang=ang+wkc(j,i)
  enddo
enddo
eke=f12*eke*dsumi
ens=f12*ens*dsumi
ang=    ang*dsumi

 !Record diagnostics in monitor.asc:
write(20,'(f12.5,3(1x,f14.9))') t,eke,ens,ang
write(* ,'(f12.5,3(1x,f14.9))') t,eke,ens,ang

!------------------------------------------------------------
 !Increment counter for dumping vorticity field:
idump=idump+1

 !Define single precision values for output:
do i=1,nt
  do j=1,ng
    zzr4(j,i)=real(zz(j,i))
  enddo
enddo
tr4=real(t)

 !Write field for later imaging or diagnostics:
write(40,rec=idump) tr4,zzr4

!------------------------------------------------------------------------
 !Write contours to the cont subdirectory:
if (loop .gt. 0) then
  ofile='000'
  write(ofile(1:3),'(i3.3)') loop

  write(80,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt

  open(81,file='cont/index'//ofile,form='unformatted',status='replace')
  write(81) np(1:n),i1(1:n),ind(1:n)
  close(81)

  open(82,file='cont/nodes'//ofile,form='unformatted',status='replace')
  write(82) x(1:npt),y(1:npt),z(1:npt)
  close(82)
endif

return
end subroutine

!=======================================================================

subroutine finalise

implicit none

 !Close output files (opened in subroutine initialise):
close(20)
close(30)
close(40)
close(80)

write(*,*) 'casl completed normally' 

return 
end subroutine

 !End main program
end program
!=======================================================================
