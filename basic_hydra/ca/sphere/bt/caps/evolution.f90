module evolution

! Module contains subroutines to evolve fields as detailed in spe.f90.

use common

implicit none

 !Full PV anomaly (relative vorticity) & part associated with the contours:
double precision:: qq(ng,nt),qc(ng,nt)

 !Integrating factors used for residual PV evolution:
double precision:: eph(ng,nt),emh(ng,nt),epf(ng,nt),emf(ng,nt)

 !Velocity field and contour velocities:
double precision:: uu(0:ng+1,nt),vv(0:ng+1,nt)
double precision:: u(npm),v(npm),w(npm)

 !Variables
double precision:: twist

contains 

!=============================================================
subroutine advect

! Main subroutine for advecting fields and contours

implicit none

 !Local variables:
double precision,parameter:: twistmax=2.5
 !twistmax: the maximum value of the time integral of |zeta|_max
 !          between regularisations of the contours.
integer,parameter:: nregmax=20
 !nregmax: every nregmax regularisations, the code stops to rebuild
 !         the contours in a separate module (congen.f90).
integer:: nreg

!--------------------------------------------------------------
 !Initialise residual PV field after recontouring:
call init

 !Initialise counter for regularisations:
nreg=0

 !Used for regularising contours:
twist=zero

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (nint(tmax/tsim) .le. nperiod) 
  do while (t .lt. tmax) 

     !Evolve contours and fields from t to t + dt:
    call advance

    if (twist .gt. twistmax) then
       !Check whether to regularise contours or rebuild them
      nreg=nreg+1

      if (nreg .eq. nregmax) then
         !Rebuild the contours using module congen:
        call prepare
        return
      endif
 
       !Regularise the contours (surgery + node redistribution):
      call surgery
      twist=twist-twistmax
    endif
  enddo  

   !Write data for later post-processing:
  call writedata(1)

   !Continue for another simulation period:
  tmax=t+tsim
enddo
 !End time loop 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !Reaching this point, the code has completed the simulation

 !If t is a multiple of tsim, write data for later post-processing:
if (abs(t/tsim-dble(nint(t/tsim))) .lt. 1.d-6) call writedata(1)

 !Set t to final time to ensure main code loop stops:
t=tfin

return
end subroutine

!=======================================================================

subroutine init
! Initialises PV fields after recontouring and before timestepping

implicit double precision(a-h,o-z)
implicit integer(i-n)

!------------------------------------------------------------------------
 !Update gridded version of contour PV (absolute vorticity):
call con2grid(qc)

 !Define relative vorticity associated with the contours:
do i=1,nt
  do j=1,ng
    qq(j,i)=qc(j,i)-cof(j)
  enddo
enddo

 !FFT qq in longitude (semi-spectral):
call forfft(ng,nt,qq,trig,factors) 
 !qq contains the PV anomaly arising from the contours in 
 !semi-spectral space (called "qc" in the comments below).
 !(qs contains the entire PV anomaly in semi-spectral space)

 !Calculate the residual PV from qd = (1-F)*(q-qc):
do m=1,nt
  do j=1,ng
    qd(j,m)=qs(j,m)-qq(j,m)
    qq(j,m)=qs(j,m)
  enddo
enddo

 !Apply high-pass Butterworth filter 1-F to (q-qc) in array qd:
call hifilter(qd)

!--------------------------------------------------------------------
 !Time which we are going to compute until:
tmax=tsim*(one+dble(int(t/tsim)))
 !Here, t is the initial time, read in above.

dt=dtmax
 !The time step is reset in subroutine adapt but needs to be defined

 !Write initial conditions for later post-processing:
if (t .eq. zero) then
  call writedata(0)
endif

return
end subroutine

!=======================================================================

subroutine prepare
! Prepares PV fields for recontouring

implicit double precision(a-h,o-z)
implicit integer(i-n)

!-------------------------------------------
 !Re-initialise qs & qd before recontouring:
 !    Reset qs = q = F*qs + (1-F)*qc + qd
 !      and qd = (1-F)*(q-qc)
call reset
 !qq contains the semi-spectral PV anomaly after this call.

 !Get qq in physical space:
call revfft(ng,nt,qq,trig,factors)

 !Define vorticity for recontouring:
do i=1,nt
  do j=1,ng
    qd(j,i)=qq(j,i)-qc(j,i)+cof(j)
  enddo
enddo
 !Note: qc is the PV associated with the contours.

return
end subroutine

!=======================================================================

subroutine advance
! Integrates the equations of motion from time t to time t + dt,
! where dt is dynamically adapted in subroutine adapt.
! *** Uses the 4th-order Runge-Kutta method ***

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local variables:
double precision:: xi(npm),yi(npm),zi(npm)
double precision:: xf(npm),yf(npm),zf(npm)
double precision:: qsi(ng,nt),qsf(ng,nt)
double precision:: qdi(ng,nt),qdf(ng,nt)
double precision:: sqs(ng,nt),sqd(ng,nt)

!------------------------------------------------------------------------
 !Re-initialise qs & qd at the beginning of the time step:
 !    Reset qs = q = F*qs + (1-F)*qc + qd
 !      and qd = (1-F)*(q-qc)
call reset
 !qq contains the semi-spectral PV anomaly after this call.
 !qd contains the semi-spectral residual PV.

!-----------------------------------------------------------------------
 !Find gridded velocity, possibly adapt the time step, interpolate the
 !velocity at each contour node, and compute the qd & qd source terms:
call invert(sqs,sqd,0)

!---------------------------------------------------------------------
 !RK4 predictor step to time t0 + dt/2:

 !Contour nodes:
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

 !Semi-spectral fields:
do m=1,nt
  do j=1,ng
    qsi(j,m)=qs(j,m)
    qs(j,m)=emh(j,m)*(qsi(j,m)+dt2*sqs(j,m))
    qsf(j,m)=qsi(j,m)+dt6*sqs(j,m)
    qdi(j,m)=qd(j,m)
    qd(j,m)=emh(j,m)*(qdi(j,m)+dt2*sqd(j,m))
    qdf(j,m)=qdi(j,m)+dt6*sqd(j,m)
  enddo
enddo

!---------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt/2:
t=t+dt2
call invert(sqs,sqd,1)

 !Contour nodes:
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

 !Semi-spectral fields:
do m=1,nt
  do j=1,ng
    qs(j,m)=emh(j,m)*(qsi(j,m)+dt2*sqs(j,m))
    qsf(j,m)=qsf(j,m)+dt3*sqs(j,m)
    qd(j,m)=emh(j,m)*(qdi(j,m)+dt2*sqd(j,m))
    qdf(j,m)=qdf(j,m)+dt3*sqd(j,m)
  enddo
enddo

!---------------------------------------------------------------------
 !RK4 predictor step at time t0 + dt:
call invert(sqs,sqd,1)

 !Contour nodes:
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

 !Semi-spectral fields:
do m=1,nt
  do j=1,ng
    qs(j,m)=emf(j,m)*(qsi(j,m)+dt*sqs(j,m))
    qsf(j,m)=qsf(j,m)+dt3*sqs(j,m)
    qd(j,m)=emf(j,m)*(qdi(j,m)+dt*sqd(j,m))
    qdf(j,m)=qdf(j,m)+dt3*sqd(j,m)
  enddo
enddo

!---------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt:
t=t+dt2
call invert(sqs,sqd,2)

 !Contour nodes:
do i=1,npt
  xnew=xf(i)+dt6*u(i)
  ynew=yf(i)+dt6*v(i)
  znew=zf(i)+dt6*w(i)
  fac=one/sqrt(xnew**2+ynew**2+znew**2)
  x(i)=fac*xnew
  y(i)=fac*ynew
  z(i)=fac*znew
enddo

 !Semi-spectral fields:
do m=1,nt
  do j=1,ng
    qs(j,m)=emf(j,m)*(qsf(j,m)+dt6*sqs(j,m))
    qd(j,m)=emf(j,m)*(qdf(j,m)+dt6*sqd(j,m))
  enddo
enddo

 !Apply latitudinal hyperviscous damping:
call latdamp(qd)
call latdamp(qs)

!-----------------------------------------------------------------------
 !Add any stochastic forcing to qd here:
if (stoch) then
   !rms value of added vorticity:
  dzrms=sqrt(two*esr*dt)
   !generate random field with this rms value and add to qd:
  call ranspec(qd,dzrms)
endif

return
end subroutine

!==========================================================================

subroutine invert(sqs,sqd,lev)
! Calculates the velocity field (uu,vv) given the PV anomaly field (qq).  
! The velocity field is then interpolated at all nodes as (u,v,w).
! The nonlinear source terms (sqs & sqd) for qs & qd are also calculated.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: sqs(ng,nt),sqd(ng,nt)

 !Work arrays:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt)
double precision::  zz(ng,nt)

!------------------------------------------------------------------------
 !Get the semi-spectral PV anomaly (qq):
if (lev .gt. 0) call genpv
 !Note: if lev = 0, qq has already been obtained in subroutine reset.

 !Invert Laplace's operator on the PV anomaly:
call laplinv(qq,wka,wkb)
 !Here the streamfunction psi is wka while wkb = tau*d(psi)/dlat.

 !Compute d(psi)/dlon = wkc:
call deriv(ng,nt,rk,wka,wkc)

 !Get physical space velocity:
call revfft(ng,nt,wkb,trig,factors)
call revfft(ng,nt,wkc,trig,factors)  

 !Copy into uu & vv which include polar points (see velint in contours.f90):
 !and correct meridional velocity by dividing by cos(lat):
do i=1,nt
  do j=1,ng
    uu(j,i)=      -wkb(j,i)
    vv(j,i)=tdr(j)*wkc(j,i)
  enddo
enddo
 !Above, tdr = tau/rho.

 !Interpolate the velocity at the contour nodes:
call velint(uu,vv,u,v,w)

!----------------------------------------------------------
 !Save various fields approximately every tsave time units:
if (lev .eq. 0) then

   !Get qq in physical space (as zz, the relative vorticity):
  do m=1,nt
    do j=1,ng
      zz(j,m)=qq(j,m)
    enddo
  enddo
  call revfft(ng,nt,zz,trig,factors)

  itime=int(t/tsave+small)
  tcrit=tsave*dble(itime)
  if ((t-dt-tcrit)*(t+small*tsave-tcrit) .lt. zero) then
    if (itime .gt. iref) then
      call dump(zz)
      iref=itime
    endif
  endif

   !Possibly adapt the time step and adjust damping rates:
  call adapt(zz)
endif

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !Get source term (sqs) for spectral PV field, qs_t:
do m=1,nt
  do j=1,ng
    wkc(j,m)=qs(j,m)
  enddo
enddo

 !Calculate longitudinal derivative spectrally:
call deriv(ng,nt,rk,wkc,wka)
 !Return to physical space:
call revfft(ng,nt,wka,trig,factors) 
 !Calculate latitudinal derivative spectrally:
call latder(wkc,wkb)
 !wkb is returned in physical space here

 !Get NL tendency for qs:
do i=1,nt
  do j=1,ng
    sqs(j,i)=-uu(j,i)*rhoi(j)*wka(j,i)-vv(j,i)*(wkb(j,i)+bet(j))
  enddo
enddo

 !FFT sqs in longitude after de-aliasing:
call dealiase(sqs)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !Get source term (sqd) for spectral residual PV field, qd_t:
do m=1,nt
  do j=1,ng
    wkc(j,m)=qd(j,m)
  enddo
enddo

 !Calculate longitudinal derivative spectrally:
call deriv(ng,nt,rk,wkc,wka)
 !Return to physical space:
call revfft(ng,nt,wka,trig,factors) 
 !Calculate latitudinal derivative spectrally:
call latder(wkc,wkb)
 !wkb is returned in physical space here

 !Get NL tendency for qd:
if (friction) then
   !Apply Ekman friction:
  do i=1,nt
    do j=1,ng
      sqd(j,i)=-uu(j,i)*rhoi(j)*wka(j,i)-vv(j,i)*wkb(j,i) &
             & -rekman*wke(j,i)/(one+hhp(j,i))
    enddo
  enddo
else
   !Conservative case:
  do i=1,nt
    do j=1,ng
      sqd(j,i)=-uu(j,i)*rhoi(j)*wka(j,i)-vv(j,i)*wkb(j,i)
    enddo
  enddo
endif

 !FFT sqd in longitude after de-aliasing:
call dealiase(sqd)

!-----------------------------------------------------------
if (lev .eq. 0) return

 !Apply exponential integrating factor:
if (lev .eq. 1) then
  do m=1,nt
    do j=1,ng
      sqs(j,m)=eph(j,m)*sqs(j,m)
      sqd(j,m)=eph(j,m)*sqd(j,m)
    enddo
  enddo
else
  do m=1,nt
    do j=1,ng
      sqs(j,m)=epf(j,m)*sqs(j,m)
      sqd(j,m)=epf(j,m)*sqd(j,m)
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine adapt(zz)
! Adapts the time step to ensure 
!    dt < cflmax*dl/|u|_max, 
!    dt < dtzz/|z|_max, and
!    dt < dtmax,
! where |z|_max is the maximum relative vorticity, dtzz = pi/10,
! cflmax = 0.4, dl = longitudinal grid spacing, and dtmax is 
! provided in parameters.asc.  Also adjusts the damping rate of
! residual PV.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local parameters:
double precision,parameter:: cflmax=0.7d0, dtzz=pi/10.d0
double precision,parameter:: fdtmin=0.98d0, fdtmax=one/fdtmin
 !*** The time step is adapted only if dtnew < fdtmin*dt or 
 !                                     dtnew > fdtmax*dt

 !Passed array (relative vorticity):
double precision:: zz(ng,nt)

 !Work arrays:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt)
logical change

!------------------------------------------------------------------
 !Compute energy, enstrophy and angular momentum per unit area as
 !well as max abs velocity and vorticity:
umax=small
zmax=small
do i=1,nt
  do j=1,ng
    usq=uu(j,i)**2+tsqi(j)*vv(j,i)**2
    zsq=zz(j,i)**2
    umax=max(umax,usq)
    zmax=max(zmax,zsq)
    wka(j,i)=rdt(j)*usq
    wkb(j,i)=rdt(j)*zz(j,i)**2
    wkc(j,i)=rdt(j)*rho(j)*uu(j,i)
  enddo
enddo
umax=sqrt(umax)
zmax=sqrt(zmax)

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

 !CFL time step: appears to work despite shrinking grid lengths
dtcfl=dl*cflmax/umax

 !Accurate advection time step:
dtacc=dtzz/zmax

!---------------------------------------------------------------------
 !Choose the smallest time step:
dtnew=min(dtmax,dtcfl,dtacc)

if ((dtnew .lt. fdtmin*dt) .or. (dtnew .gt. fdtmax*dt)) then
  dt=dtnew
  change=.true.
else
  change=.false.
endif

 !See how close we are to the end of the simulation:
trem=tmax-t

if (trem .lt. dt .and. trem .gt. 1.d-6*dtmax) then
  dt=trem
  change=.true.
endif

!---------------------------------------------------------------------
 !Increment the integral of max|zz|:
twist=twist+dt*zmax

 !Record diagnostics in monitor.asc:
write(18,'(f12.5,5(1x,f12.7))') t,eke,ens,ang,umax,zmax

if (.not. change) return
!----------------------------------------------------------------
 !Change the time step

 !Define Runge-Kutta time interval constants (see evolve):
dt2=f12*dt
dt3=f13*dt
dt6=f16*dt

 !Define arrays needed in spectral time stepping of qd:
zrms=sqrt(two*ens)
dfac=zrms*dt2
do m=1,nt
  do j=1,ng
    arg=min(dfac*dislon(j,m),350.d0)
    eph(j,m)=exp(arg)
    emh(j,m)=one/eph(j,m)
    epf(j,m)=eph(j,m)**2
    emf(j,m)=emh(j,m)**2
  enddo
enddo

 !Define latitudinal hyperviscous damping:
do k=1,nt
  hyplat(k)=one/(one+dfac*dislat(k))
enddo

return
end subroutine

!========================================================================

subroutine genpv
! Computes the PV field qq by combining the contour PV qc with that 
! in the fields qs & qd.  
! ***qq is returned in semi-spectral space***

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Work array:
double precision:: wka(ng,nt)

!-----------------------------------------------------
 !Grid the PV contours as qc:
call con2grid(qc)

 !Subtract Coriolis frequency to define anomaly (qq):
do i=1,nt
  do j=1,ng
    qq(j,i)=qc(j,i)-cof(j)
  enddo
enddo

 !FFT qq in longitude:
call forfft(ng,nt,qq,trig,factors) 
 !qq contains the PV anomaly arising from the contours in 
 !semi-spectral space (called "qc" in the comments below)

 !Store qs-qc in wka in preparation for filtering:
do m=1,nt
  do j=1,ng
    wka(j,m)=qs(j,m)-qq(j,m)
  enddo
enddo

 !Apply low-pass filter F to qs-qc:
call lofilter(wka)
 !wka = F(qs-qc) (now semi-spectral)

 !Define PV anomaly from qq = F(qs-qc)+qd+qc:
do m=1,nt
  do j=1,ng
    qq(j,m)=wka(j,m)+qd(j,m)+qq(j,m)
  enddo
enddo

 !Remove global mean value of qq (now the full PV anomaly)
avqq=f1112*(qq(1,1)*rdt(1)+qq(ng,1)*rdt(ng))
do j=2,ngm1
  avqq=avqq+qq(j,1)*rdt(j)
enddo
avqq=avqq*rsumi
do j=1,ng
  qq(j,1)=qq(j,1)-avqq
enddo

return
end subroutine

!========================================================================

subroutine reset
! Resets qs & qd at the beginning of each time step, or
! just before recontouring at the end of a run.  Also 
! defines the full PV anomaly qq, in semi-spectral space.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Work array:
double precision:: wka(ng,nt)

!------------------------------------------------------------------------
 !Grid the PV contours as qc:
call con2grid(qc)

 !Subtract Coriolis frequency to define anomaly (qq):
do i=1,nt
  do j=1,ng
    qq(j,i)=qc(j,i)-cof(j)
  enddo
enddo

 !FFT qq in longitude:
call forfft(ng,nt,qq,trig,factors) 
 !qq contains the PV anomaly arising from the contours in 
 !semi-spectral space (called "qc" in the comments below)

 !Store qs-qc in wka in preparation for filtering:
do m=1,nt
  do j=1,ng
    wka(j,m)=qs(j,m)-qq(j,m)
  enddo
enddo

 !Apply low-pass filter F to qs-qc:
call lofilter(wka)
 !wka = F(qs-qc) (now semi-spectral)

 !Add on qd and obtain reset value of qs <-- F(qs-qc)+qd+qc:
do m=1,nt
  do j=1,ng
    qd(j,m)=wka(j,m)+qd(j,m)
    qs(j,m)= qd(j,m)+qq(j,m)
    qq(j,m)=         qs(j,m)
  enddo
enddo
 !Note: now qd = *new* qs-qc

 !Reset qd = (1-F)(qs-qc):
call hifilter(qd)

 !Remove global mean value of qq (now the full PV anomaly)
avqq=f1112*(qq(1,1)*rdt(1)+qq(ng,1)*rdt(ng))
do j=2,ngm1
  avqq=avqq+qq(j,1)*rdt(j)
enddo
avqq=avqq*rsumi
do j=1,ng
  qq(j,1)=qq(j,1)-avqq
enddo
 !qq is in semi-spectral space.

return
end subroutine

!========================================================================

subroutine dump(zz)
! Writes the relative vorticity field (zz) to the main job directory

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed array:
double precision:: zz(ng,nt)

 !Local variables & work arrays:
real:: zzr4(ng,nt),tr4

!------------------------------------------------------------
 !Increment counter for dumping data:
idump=idump+1

 !Define single precision values of various fields for output:
do i=1,nt
  do j=1,ng
    zzr4(j,i)=real(zz(j,i))
  enddo
enddo
tr4=real(t)

 !Write fields for later imaging or diagnostics:
write(45,rec=idump) tr4,zzr4

return
end subroutine

!========================================================================

subroutine writedata(iopt)
! This routine writes the current contours and the residual 
! PV to the cont subdirectory and the relative vorticity zz 
! to the grid subdirectory.

! If iopt = 1, qd is recomputed by a call to reset first (this
! is not necessary at t = 0 in subroutine init).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local variables:
double precision:: zz(ng,nt)
real:: zzr4(ng,nt),tr4
character(len=3):: ofile

!----------------------------------------------------------------
 !Reset qs = q = F*qs + (1-F)*qc + qd  &  qd = (1-F)*(q-qc):
if (iopt .eq. 1) call reset
 !qq contains the full gridded PV anomaly after this call.
 !qd contains the *spectral* residual PV.

 !Get qq in physical space (as zz, the relative vorticity):
do m=1,nt
  do j=1,ng
    zz(j,m)=qq(j,m)
  enddo
enddo
call revfft(ng,nt,zz,trig,factors)

 !Define residual needed for constructing ultra-fine PV field
 !in post processing:
do i=1,nt
  do j=1,ng
    zzr4(j,i)=real(zz(j,i)-qc(j,i)+cof(j))
  enddo
enddo
 !Note: qc is the total PV coming from the contours.

!------------------------------------------------------------------------
 !Write data:
loop=nint(t/tsim)
irec=loop+1
ofile='000'
write(ofile(1:3),'(i3.3)') loop

tr4=real(t)
 !Write contours to the cont subdirectory:
write(80,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt

open(81,file='cont/index'//ofile,form='unformatted',status='replace')
write(81) np(1:n),i1(1:n),ind(1:n)
close(81)

open(82,file='cont/nodes'//ofile,form='unformatted',status='replace')
write(82) x(1:npt),y(1:npt),z(1:npt)
close(82)

 !Write residual PV contained in zzr4 to cont subdirectory:
write(83,rec=irec) tr4,zzr4

 !Write relative vorticity zz to grid subdirectory
 !Define single precision values for output:
do i=1,nt
  do j=1,ng
    zzr4(j,i)=real(zz(j,i))
  enddo
enddo

 !Write fields for later imaging or diagnostics:
write(35,rec=irec) tr4,zzr4

 !Write final data for imaging:
if (loop .eq. nperiod) call dump(zz)

return
end subroutine

!=======================================================================

 !Main end module
end module
