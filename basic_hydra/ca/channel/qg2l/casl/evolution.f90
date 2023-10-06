module evolution

! Module contains subroutines to evolve PV according to the 
! algorithm detailed in casl.f90.

use common

implicit none

 !Energies:
double precision:: ekepre,ekepost,apepre,apepost

 !Physical fields:
double precision::    qq(0:ny,0:nxm1,nz),qqpre(0:ny,0:nxm1,nz)
double precision::   dd1(0:ny,0:nxm1),  dd1pre(0:ny,0:nxm1)
double precision:: qspre(0:ny,0:nxm1,nz),qdpre(0:ny,0:nxm1,nz)
double precision::    uu(0:ny,0:nxm1,nz),   vv(0:ny,0:nxm1,nz)
double precision::    pp(0:ny,0:nxm1,nz),   qc(0:ny,0:nxm1,nz)

 !Semi-Lagrangian advection arrays:
double precision::  x0(0:ny,0:nxm1,nz), y0(0:ny,0:nxm1,nz)
double precision:: ula(0:ny,0:nxm1,nz),vla(0:ny,0:nxm1,nz)

!Internal subroutine definitions (inherit global variables):

contains 

!=============================================================
subroutine advect

! Main subroutine for advecting fields and contours

implicit none

 !Local variables:
double precision,parameter:: twistmax=2.5d0
!      twistmax: the maximum value of the time integral of |zeta|_max
!                between regularisations of the contours.
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops 
!      to rebuild the contours in a separate memory space.
integer:: ireg,igsave,icsave,ix,iy,iz

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Used for regularising contours:
twist=zero

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Perform surgery & field reset 
  if (twist .gt. twistmax) then
    ireg=ireg+1
     !Don't continue if maximum number of regularisations reached;
     !it is time to recontour (pass control back to main program):
    if (ireg .eq. nregmax) then
       !Prepare for recontouring:
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise contours (surgery + node redistribution):
    call surgery
     !Record contour complexity to complexity.asc:
    write(14,'(1x,f13.5,1x,i8,1x,i9)') t,nq,nptq

     !Convert PV contours to gridded values (qc):
    call con2grid(qc)
     !Reset qs and qd:
    call reset(qc,qs,qd,qavg)

     !Copy gridded PV fields to old time level:
    call copyfields

     !Update twist parameter:
    twist=twist-twistmax
  endif

   !Adjust timestep (dt) on maximum vorticity magnitude:
  call adapt(igsave,icsave)
   !Advect PV from time t to t + dt:
  call advance

   !Update the time:
  t=t+dt

   !Possibly save PV & energy at chosen save time (tgrid):
  if (igsave .eq. 1) call savegrid

   !Possibly save contours and residual PV (qd) for post processing:
  if (icsave .eq. 1) call savecont
  
   !Copy new fields into previous time:
  call copyfields
enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
return
end subroutine

!=======================================================================

subroutine prepare

! Prepares fields and contour intervals for recontouring

implicit none

double precision:: qqmin(nz),qqmax(nz)
integer:: ix,iy,iz

 !Prepare PV residual for recontouring:
do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      qd(iy,ix,iz)=qq(iy,ix,iz)-qc(iy,ix,iz)
      qs(iy,ix,iz)=qq(iy,ix,iz)
    enddo
  enddo
enddo

 !Update also PV jumps:
call contint(qs,ncontq,qjump,qqmin,qqmax)

 !Write information to log file:
do iz=1,nz
  write(*,'(a,i1,a,1x,f13.5,1x,f11.5,1x,f9.5)') ' layer',iz, &
              & ' q_min, q_max, qjump = ',qqmin(iz),qqmax(iz),qjump(iz)
enddo
 
return
end subroutine

!=======================================================================

subroutine copyfields

! Copies new fields into previous time

implicit none

integer:: ix,iy,iz

if (damping) then
   !Copy PV source term:
  do iz=1,nz
    do ix=0,nxm1
      do iy=0,ny
        sqpre(iy,ix,iz)=sq(iy,ix,iz)
      enddo
    enddo
  enddo
endif

 !Copy PV fields:
do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      qqpre(iy,ix,iz)=qq(iy,ix,iz)
      qdpre(iy,ix,iz)=qd(iy,ix,iz)
      qspre(iy,ix,iz)=qs(iy,ix,iz)
    enddo
  enddo
enddo

 !Copy interface displacements (multiplied by f_0/(H1+H2)):
if (barot) then
   !Here dd2 = 0
  do ix=0,nxm1
    do iy=0,ny
      dd1pre(iy,ix)=dd1(iy,ix) 
    enddo
  enddo
else
  do ix=0,nxm1
    do iy=0,ny
      dd1pre(iy,ix)=dd1(iy,ix) 
      dd2pre(iy,ix)=dd2(iy,ix) 
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine init

! Initialises quantities needed for normal time integration following 
! contour regeneration

implicit none

integer:: ix,iy,iz

!---------------------------------------------
 !Record contour complexity to complexity.asc:
write(14,'(1x,f13.5,1x,i8,1x,i9)') t,nq,nptq

!---------------------------------------------------
 !Convert PV contours to gridded values (qc):
call con2grid(qc)

 !Define residual PV qd = qs-qc-F[qs-qc]
do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      ula(iy,ix,iz)=qs(iy,ix,iz)-qc(iy,ix,iz)
      vla(iy,ix,iz)=ula(iy,ix,iz)
    enddo
  enddo
enddo

call filter(ula,0,2)

 !Finish qd calculation & copy gridded PV to old time level:
do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      qd(iy,ix,iz)=vla(iy,ix,iz)-ula(iy,ix,iz)
      qqpre(iy,ix,iz)=qs(iy,ix,iz)
      qdpre(iy,ix,iz)=qd(iy,ix,iz)
      qspre(iy,ix,iz)=qs(iy,ix,iz)
    enddo
  enddo
enddo

!--------------------------------------------------------
 !Get the initial velocity (uu,vv) and streamfunction pp:
call inversion

 !Calculate the relaxation terms for PV (sqpre) if thermal/Ekman 
 !damping is active:
if (damping) call relax(sqpre)
 !Note: qqpre and sqpre are needed by subroutine advance.

 !Copy interface displacements (multiplied by f_0/(H1+H2)):
if (barot) then
   !Here dd2 = 0
  do ix=0,nxm1
    do iy=0,ny
      dd1pre(iy,ix)=dd1(iy,ix) 
    enddo
  enddo
else
  do ix=0,nxm1
    do iy=0,ny
      dd1pre(iy,ix)=dd1(iy,ix) 
      dd2pre(iy,ix)=dd2(iy,ix) 
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine relax(qqsrc)

! Computes the thermal and/or Ekman damping terms needed in the PV 
! evolution equation.

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! *** If any source term were to change the mean PV, then qavg would 
! *** need to be updated here.
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed array:
double precision:: qqsrc(0:ny,0:nxm1,nz)
 !Local array:
double precision:: zz1(0:ny,0:nxm1)

!-----------------------------------------------------------------
if (heating) then
   !Apply thermal damping:
  if (barot) then
     !Upper interface displacement is zero here
    do ix=0,nxm1
      do iy=0,ny
        difdd1=dd1(iy,ix)-dd1eq(iy,ix)
        qqsrc(iy,ix,1)= therm1*difdd1
        qqsrc(iy,ix,2)=-therm2*difdd1
      enddo
    enddo
  else
    do ix=0,nxm1
      do iy=0,ny
        difdd1=dd1(iy,ix)-dd1eq(iy,ix)
        difdd2=dd2(iy,ix)-dd2eq(iy,ix)
        qqsrc(iy,ix,1)=therm1*difdd1
        qqsrc(iy,ix,2)=therm2*(difdd2-difdd1)
      enddo
    enddo
  endif
else
   !Initialise qqsrc:
  do ix=0,nxm1
    do iy=0,ny
      qqsrc(iy,ix,1)=zero
      qqsrc(iy,ix,2)=zero
    enddo
  enddo
endif

if (friction) then 
   !Apply Ekman damping, -rekman*zeta_1, to lowest layer:
  if (barot) then
    do ix=0,nxm1
      do iy=0,ny
        zz1(iy,ix)=qq(iy,ix,1)-bety(iy)-fhb(iy,ix)+h1inv*dd1(iy,ix)
      enddo
    enddo
  else
    do ix=0,nxm1
      do iy=0,ny
        zz1(iy,ix)=qq(iy,ix,1)-bety(iy)-fhb(iy,ix)+h1inv*dd1(iy,ix)
      enddo
    enddo
  endif
  do ix=0,nxm1
    do iy=1,nym1
      qqsrc(iy,ix,1)=qqsrc(iy,ix,1)-rekman*zz1(iy,ix)
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine pvforce

! Adds stochastic forcing (in the form of hetons) to the PV (in qd)
! at the end of a time step.  The forcing is purely applied to mode 2
! (the baroclinic mode with the largest kd).  A compensating uniform
! PV is added to ensure no net PV is input into either layer - this is 
! needed to conserve the mass of each layer.

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! *** If the forcing were to change the mean PV, then qavg would need
! *** to be updated here.
! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

implicit none

 !Local variables:
double precision:: dqf2(0:nyfp1,0:nxfm1),dq2(0:ny,0:nxm1)
double precision:: dnx,dny,xx,xv,yy,yv,delx,dely,avgdq2,dq2p
integer:: nvor,k,ix,iy,ix1,ix2,iy1,iy2,mix

!-----------------------------------------------------------------
 !Add circular vortex hetons here to force the flow (they are 
 !converted to gridded PV values and added to qd):

 !Add nvor vortices (as +/- pairs in lower/upper layers):
nvor=2*nint(f12*(dnvor*t-totnvor))
if (nvor .eq. 0) return

 !Initialise array to uptake added PV in vertical mode 2:
do ix=0,nxfm1
  do iy=0,nyf
    dqf2(iy,ix)=zero
  enddo
enddo

 !The vortices are place randomly in the domain:
dnx=dble(nxf)-small
dny=dble(nyf)-small
do k=1,nvor
   !Place vortices at random in the domain and work out which 
   !grid box centres lie within the vortices:
  xx=dnx*rand(0)
  xv=glxf*xx
  ix1=nint(nxf+xx-rhetgx)-nxf
  ix2=nint(xx+rhetgx)

  yy=dny*rand(0)
  yv=glyf*yy
  iy1=max(0,nint(nyf+yy-rhetgy)-nyf)
  iy2=min(nyf,nint(yy+rhetgy))

  do ix=ix1,ix2
    mix=mod(nxf+ix,nxf)
    delx=xv-glxf*(dble(ix)+f12)
    do iy=iy1,iy2
      dely=yv-glyf*(dble(iy)+f12)
      if (delx**2+dely**2 .lt. rhetsq) dqf2(iy,mix)=dqf2(iy,mix)+vorvor
    enddo
  enddo
enddo

 !Coarsen to inversion grid:
call coarsen(dqf2,dq2)

 !Total number of vortices/dipoles added so far:
totnvor=totnvor+dble(nvor)

 !Remove average dq2:
call average(dq2,avgdq2)

 !Add dq2 to qd so that the upper layer is anti-cyclonic and the lower
 !layer is cyclonic:
do ix=0,nxm1
  do iy=0,ny
    dq2p=dq2(iy,ix)-avgdq2
    qd(iy,ix,1)=qd(iy,ix,1)-vect12*dq2p
    qd(iy,ix,2)=qd(iy,ix,2)-vect22*dq2p
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine inversion

! Combines the contour and grid PV fields then inverts the PV to
! obtain the velocity and streamfunction.  Also calculates the
! (scaled) interface displacements, dd1 & dd2

implicit none
integer:: ix,iy

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qc)

 !Combine fields to update qq with full field:
call combine(qq,qc,qs,qd,qavg)

 !Invert PV to obtain velocity field (uu,vv) and streamfunction (pp):
call main_invert(qq,fhb,t,uu,vv,pp)
 !Note: fhb contains f_0*H_b/H_1 where H_b is the bottom topography

 !Define the scaled interface displacements, f_0*eta_j/(H_1+H_2):
if (barot) then
   !Upper interface is flat in this case (dd2 = 0):
  do ix=0,nxm1
    do iy=0,ny
      dd1(iy,ix)=h1h2kdbarsq*(pp(iy,ix,1)-pp(iy,ix,2))
    enddo
  enddo
else
  do ix=0,nxm1
    do iy=0,ny
      dd1(iy,ix)=h1h2kdbarsq*(pp(iy,ix,1)-alpha*pp(iy,ix,2))
      dd2(iy,ix)=h1h2ackdbarsq*pp(iy,ix,2)
    enddo
  enddo
endif

return
end subroutine

!=======================================================================
      
subroutine advance

! Advances PV contours and gridded fields to time t+dt by contour 
! advection and by trajectory integration using a standard 
! semi-Lagrangian (SL) scheme.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
double precision:: uq(nptq),vq(nptq),xqm(nptq),yqm(nptq)

!------------------------------------------------------------------------
 !Increments in grid units needed below for trajectory integration:
gcx=dt*glxi
gcy=dt*glyi
hgcx=f12*gcx
hgcy=f12*gcy

 !Copy current velocity field into (ula,vla) for use below then compute
 !Euler backward predictor for departure grid point (x0,y0):
do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      ula(iy,ix,iz)=uu(iy,ix,iz)
      vla(iy,ix,iz)=vv(iy,ix,iz)
      x0(iy,ix,iz)=mod(xigmax + xig(ix)-gcx*uu(iy,ix,iz),xigmax)
      y0(iy,ix,iz)=max(zero,min(yig(iy)-gcy*vv(iy,ix,iz),yigmax))
    enddo
  enddo
enddo
 !Note, (uu,vv) is used since we have no other velocity field available

 !Prepare contour evolution; get velocity on PV contour nodes:
if (nptq .gt. 0) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xq(i)+hfdt*uq(i)
    xqm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqm(i)=min(ymax,max(ymin,yq(i)+hfdt*vq(i)))
    xx=xq(i)+dt*uq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=min(ymax,max(ymin,yq(i)+dt*vq(i)))
  enddo
endif

 !Iterate to converge on implicit trapezoidal integration:
do iter=1,niter
   !Obtain qs & qd at time t + dt:
  call sl_step

   !Obtain qq & hence uu & vv at time t + dt:
  call inversion

   !Correct departure grid point (x0,y0):
  do iz=1,nz
    do ix=0,nxm1
      do iy=0,ny
         !Obtain old velocity (time t) at the departure grid point using
         !bi-linear interpolation of (ula,vla):
        ix0=int(x0(iy,ix,iz))
        ix1=ixp(ix0)
        px=x0(iy,ix,iz)-dble(ix0)
        pxc=one-px

        iy0=int(y0(iy,ix,iz))
        iy1=iyp(iy0)
        py=y0(iy,ix,iz)-dble(iy0)
        pyc=one-py

        uod=pyc*(pxc*ula(iy0,ix0,iz)+px*ula(iy0,ix1,iz)) &
      &     +py*(pxc*ula(iy1,ix0,iz)+px*ula(iy1,ix1,iz))

        vod=pyc*(pxc*vla(iy0,ix0,iz)+px*vla(iy0,ix1,iz)) &
      &     +py*(pxc*vla(iy1,ix0,iz)+px*vla(iy1,ix1,iz))

        x0(iy,ix,iz)=mod(xigmax + xig(ix)-hgcx*(uod+uu(iy,ix,iz)),xigmax)
        y0(iy,ix,iz)=max(zero,min(yig(iy)-hgcy*(vod+vv(iy,ix,iz)),yigmax))
         !(uu,vv) is the new velocity (time t+dt) at the arrival grid point
      enddo
    enddo
  enddo

   !Update the PV contour points:
  if (nptq .gt. 0) then
    call velint(uu,vv,uq,vq)
    do i=1,nptq
      xx=xqm(i)+hfdt*uq(i)
      xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
      yq(i)=min(ymax,max(ymin,yqm(i)+hfdt*vq(i)))
    enddo
  endif

enddo

 !Obtain final corrected qs & qd at time t + dt:
call sl_step

 !Obtain final corrected uu & vv at time t + dt from qq:
call inversion

 !Update the PV contour points:
if (nptq .gt. 0) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqm(i)+hfdt*uq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=min(ymax,max(ymin,yqm(i)+hfdt*vq(i)))
  enddo
endif

 !Add any stochastic forcing (if present) to qd:
if (stoch) call pvforce

 !Call con2grid to get updated contour PV (qc):
call con2grid(qc)

 !Combine fields to update qq with full field:
call combine(qq,qc,qs,qd,qavg)

return
end subroutine

!=======================================================================
      
subroutine sl_step

! Interpolates qq at points (x0,y0) and adds the source integral
! from t to t+dt to qq

implicit double precision(a-h,o-z)
implicit integer(i-n)

!-------------------------------------------------------------------------
 !Integrate qs & qd using bi-cubic Lagrange interpolation of qspre & qdpre
 !at x0,y0:
call interpol(qspre,qs,x0,y0)
call interpol(qdpre,qd,x0,y0)
 !Here, we obtain only the adiabatic evolution

 !Add diabatic terms (sq) if thermal/Ekman damping is active:
if (damping) then
  call relax(sq)

  do iz=1,nz
    do ix=0,nxm1
      do iy=0,ny
         !Obtain old source (at time t) at the departure grid point using
         !bi-linear interpolation of sqpre:
        ix0=int(x0(iy,ix,iz))
        ix1=ixp(ix0)
        px=x0(iy,ix,iz)-dble(ix0)
        pxc=one-px

        iy0=int(y0(iy,ix,iz))
        iy1=iyp(iy0)
        py=y0(iy,ix,iz)-dble(iy0)
        pyc=one-py

        sqod=pyc*(pxc*sqpre(iy0,ix0,iz)+px*sqpre(iy0,ix1,iz)) &
        &    +py*(pxc*sqpre(iy1,ix0,iz)+px*sqpre(iy1,ix1,iz))

         !Integrate in time using trapezoidal rule (sq is the new source
         !at time t+dt) and add to qd:
        qd(iy,ix,iz)=qd(iy,ix,iz)+hfdt*(sqod+sq(iy,ix,iz))
      enddo
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine adapt(igsave,icsave)

! Adapts the time step to ensure dt < dtfac/max(|zeta|_max)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local parameter used for setting time step:
double precision,parameter:: dtfac=pi/40.d0
 !Local arrays for vorticity in each layer:
double precision:: zz1(0:ny,0:nxm1),zz2(0:ny,0:nxm1)
integer:: ix,iy

!------------------------------------------------------------------------
 !Compute vorticity:
if (barot) then
  do ix=0,nxm1
    do iy=0,ny
      zz1(iy,ix)=qq(iy,ix,1)-bety(iy)-fhb(iy,ix)+h1inv*dd1(iy,ix)
      zz2(iy,ix)=qq(iy,ix,2)-bety(iy)-h2inv*dd1(iy,ix)
    enddo
  enddo
else
  do ix=0,nxm1
    do iy=0,ny
      zz1(iy,ix)=qq(iy,ix,1)-bety(iy)-fhb(iy,ix)+h1inv*dd1(iy,ix)
      zz2(iy,ix)=qq(iy,ix,2)-bety(iy)+h2inv*(dd2(iy,ix)-dd1(iy,ix))
    enddo
  enddo
endif

 !Compute max abs vorticity over both layers:
zzmax=small
zzms1=zero
zzms2=zero
do ix=0,nxm1
  do iy=0,nym1
    zzmax=max(zzmax,abs(zz1(iy,ix)),abs(zz2(iy,ix)))
    zzms1=zzms1+zz1(iy,ix)**2
    zzms2=zzms2+zz2(iy,ix)**2
  enddo
enddo
zzrms1=sqrt(dsumi*zzms1)
zzrms2=sqrt(dsumi*zzms2)
 !Rms vorticity (for diagnostic purposes only):
zzrms=sqrt(dsumi*(h1m*zzms1+alpha*h2m*zzms2))
 !h1m = h1/(h1+h2*alpha) & h2m = h2*alpha/(h1+h2*alpha)

 !Time step needed for accuracy:
dtacc=dtfac/max(zzmax,srwfm)
 !The restriction on the maximum Rossby wave frequency (srwfm)
 !ensures that the fastest Rossby wave frequency is resolved.

 !Choose a new time step, limiting it to dtmax in parameters.f90:
dt=min(dtacc,dtmax)
hfdt=dt/two

 !Increment the integral of max|zz|:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
 !Record various diagnostics to monitor.asc:
write(12,'(1x,f13.5,4(1x,f14.8))') t,zzmax,zzrms1,zzrms2,zzrms

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
tgrid=tinit+tgsave*(dble(igrec)+1.d-10)
 !A small number is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  igsave=1
   !Compute kinetic & potential energy:
  call energy(ekepre,apepre)
else
   !Do not save data:
  igsave=0
endif

 !Set flag to save contour data every tcsave time units:
tcont=tinit+tcsave*dble(icrec)
if (t .le. tcont .and. t+dt .gt. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  icsave=1
else
   !Do not save data:
  icsave=0
endif

return
end subroutine

!=======================================================================

subroutine energy(eke,ape)

! Computes the kinetic and potential energy from the velocity field
! (uu,vv) and the streamfunction pp.

! Note: the energies are normalised by rho_1*(H_1+H_2), so that only
!       the density ratio alpha = rho_2/rho_1 and the depth ratios
!       h_1 = H_1/(H_1+H_2) and h_2 = H_2/(H_1+H_2).

implicit none

double precision:: eke,ape,u1l2,v1l2,u2l2,v2l2,pl1,pl2
integer:: ix,iy

 !Kinetic energy:
call l2norm(uu(0,0,1),u1l2)
call l2norm(vv(0,0,1),v1l2)
call l2norm(uu(0,0,2),u2l2)
call l2norm(vv(0,0,2),v2l2)
eke=f12*(h1*(u1l2+v1l2)+alpha*h2*(u2l2+v2l2))

 !Potential energy:
do ix=0,nxm1
  do iy=0,ny
    ula(iy,ix,1)=pp(iy,ix,1)-alpha*pp(iy,ix,2)
  enddo
enddo

call l2norm(ula(0,0,1),pl1)
call l2norm( pp(0,0,2),pl2)
ape=f12*h1h2kdbarsq*(pl1+alpha*alphac*pl2)
 !alphac = 1 - alpha here

return
end subroutine

!=======================================================================

subroutine savegrid

! Saves various quantities at the desired save time to files 
! (postprocess with r4toc2 etc...):

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: q1spec(0:max(nx,ny)),q2spec(0:max(nx,ny))
real:: q1r4(0:ny,0:nxm1),q2r4(0:ny,0:nxm1),tr4

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Needed for real*4 writes below:
tr4=real(tgrid)

 !Record to write:
irec=igrec+1

!-----------------------------------------------------
 !Store interface displacements at save time:
if (barot) then
   !Here dd2 = 0
  do ix=0,nxm1
    do iy=0,ny
      ula(iy,ix,1)=pt*dd1pre(iy,ix)+ptc*dd1(iy,ix)
    enddo
  enddo

   !Convert to real*4:
  do ix=0,nxm1
    do iy=0,ny
      q1r4(iy,ix)=real(ula(iy,ix,1))
    enddo
  enddo

  write(33,rec=irec) tr4,q1r4
else
  do ix=0,nxm1
    do iy=0,ny
      ula(iy,ix,1)=pt*dd1pre(iy,ix)+ptc*dd1(iy,ix)
      vla(iy,ix,1)=pt*dd2pre(iy,ix)+ptc*dd2(iy,ix)
    enddo
  enddo

   !Convert to real*4:
  do ix=0,nxm1
    do iy=0,ny
      q1r4(iy,ix)=real(ula(iy,ix,1))
      q2r4(iy,ix)=real(vla(iy,ix,1))
    enddo
  enddo

  write(33,rec=irec) tr4,q1r4
  write(34,rec=irec) tr4,q2r4
endif

!-----------------------------------------------------
 !Store PV in each layer at save time:
do ix=0,nxm1
  do iy=0,ny
    ula(iy,ix,1)=pt*qqpre(iy,ix,1)+ptc*qq(iy,ix,1)
    ula(iy,ix,2)=pt*qqpre(iy,ix,2)+ptc*qq(iy,ix,2)
  enddo
enddo

 !Convert to real*4:
do ix=0,nxm1
  do iy=0,ny
    q1r4(iy,ix)=real(ula(iy,ix,1))
    q2r4(iy,ix)=real(ula(iy,ix,2))
  enddo
enddo

write(31,rec=irec) tr4,q1r4
write(32,rec=irec) tr4,q2r4

!----------------------------------------------------------
 !Compute kinetic & potential energy:
call energy(ekepost,apepost)

 !Compute time interpolated energies:
apot=pt*apepre+ptc*apepost
ekin=pt*ekepre+ptc*ekepost

!-----------------------------------------------------
 !Project PV onto the two vertical modes:
do ix=0,nxm1
  do iy=0,ny
    q1p=ula(iy,ix,1)-bety(iy)
    q2p=ula(iy,ix,2)-bety(iy)
    vla(iy,ix,1)=vec11*q1p+vec12*q2p
    vla(iy,ix,2)=vec21*q1p+vec22*q2p
  enddo
enddo

 !Compute rms values (with appropriate factors which ensure that the
 !sum of the modal enstrophies equals the sum of the (mass and thickness
 !weighted) layer enstrophies):
call l2norm(vla(0,0,1),tmp)
qq1rms=sqrt(coeff1*tmp)
call l2norm(vla(0,0,2),tmp)
qq2rms=sqrt(coeff2*tmp)

 !Write diagnostics to the files norms.asc & ene.asc:
write(13,'(f8.2,2(1x,f14.9))') tgrid,qq1rms,qq2rms
write(15,'(f8.2,3(1x,f14.9))') tgrid,ekin,apot,ekin+apot

write(*,'(a,f8.2,3(a,f11.7))') &
    & ' t = ',tgrid,'  K = ',ekin,'  P = ',apot,'  E = K+P = ',ekin+apot

 !Compute the 1d PV spectrum for each vertical mode:
call spec1d_fc(vla(0,0,1),q1spec)
call spec1d_fc(vla(0,0,2),q2spec)

sq1spec=zero
sq2spec=zero
do k=0,kmax
  q1spec(k)=coeff1*q1spec(k)
  q2spec(k)=coeff2*q2spec(k)
  sq1spec=sq1spec+q1spec(k)
  sq2spec=sq2spec+q2spec(k)
   !Normalise to take into account uneven sampling of wavenumbers
   !in each shell [k-1/2,k+1/2]:
  q1spec(k)=spmf(k)*q1spec(k)
  q2spec(k)=spmf(k)*q2spec(k)
enddo
sq1spec=8.d0*sq1spec*dsumi
sq2spec=8.d0*sq2spec*dsumi
 !Write out spectrum to file:
write(51,'(f8.2,4(1x,e14.7),1x,i5)') tgrid,sq1spec,qq1rms**2, & 
                             &             sq2spec,qq2rms**2,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(51,'(3(1x,f12.8))') alk(k),log10(q1spec(k)),log10(q2spec(k))
enddo
 !Note: alk(k) = log_10(k)

!--------------------------------------------------------------
 !Increment counter for direct file access:
igrec=igrec+1

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves PV contours and residual PV for use by congen.f90 & diagnostics

implicit double precision(a-h,o-z)
implicit integer(i-n)

real:: qdr4(0:ny,0:nxm1),tr4
integer:: iop(nq)
character(len=3):: pind

write(*,'(a,f13.5)') ' Saving contours at t = ',t

write(pind(1:3),'(i3.3)') icrec

tr4=real(t)

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f13.5,4(1x,f16.12))') nq,nptq,t, &
                    & qjump(1),qjump(2),qavg(1),qavg(2)

 !Save PV contours if any exist:
if (nq .gt. 0) then
   !First form iop; open/closed indicator:
  do j=1,nq
    iop(j)=nextq(i2q(j))/i1q(j)
     !iop = 0 for an open contour, and 1 for a closed one
  enddo
  open(81,file='contours/qqindex'//pind,form='unformatted', &
      & access='direct',status='replace',recl=20*nq)
  write(81,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq),layq(1:nq),iop(1:nq)
  close(81)

  open(82,file='contours/qqnodes'//pind,form='unformatted', &
      & access='direct',status='replace',recl=16*nptq)
  write(82,rec=1) xq(1:nptq),yq(1:nptq)
  close(82)
endif

 !Save residual needed to build ultra-fine-grid PV with congen:
open(83,file='contours/qqresi'//pind,form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
do iz=1,nz
  do ix=0,nxm1
    do iy=0,ny
      qdr4(iy,ix)=real(qq(iy,ix,iz)-qc(iy,ix,iz))
    enddo
  enddo
  write(83,rec=iz) tr4,qdr4
enddo
close(83)

!--------------------------------------------------------------
 !Increment counter for naming next direct-access output files:
icrec=icrec+1

return
end subroutine

!=======================================================================

 !Main end module
end module
