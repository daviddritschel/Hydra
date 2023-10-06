module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

double precision:: qt(ng,ng),qc(ng,ng),qd(ng,ng) !Various PVs
double precision:: dd(ng,ng) !Physical space velocity divergence

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
integer:: ireg,itime,jtime
logical:: ggen

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

   !Save data periodically:
  itime=nint(t/dt)
  jtime=itime/ngsave
  if (ngsave*jtime .eq. itime) then
    call inversion
    call savegrid(jtime+1)
    ggen=.false.
  else
    ggen=.true.
  endif
   !ggen is used to indicate if calling inversion is needed in advance below
  jtime=itime/ncsave
  if (ncsave*jtime .eq. itime) call savecont(jtime+1)

   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

     !Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax) then
       !Prepare PV residual qr for recontouring (and preserve qs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record active contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i9,1x,i10)') t,nq,nptq

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance(ggen)
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Possibly save final data:
itime=nint(t/dt)
jtime=itime/ngsave
if (ngsave*jtime .eq. itime) then
  call inversion
  call savegrid(jtime+1)
endif
jtime=itime/ncsave
if (ncsave*jtime .eq. itime) call savecont(jtime+1)

return
end subroutine

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration

implicit none

 !Local variables:
double precision:: qa(ng,ng)

!------------------------------------------------------------------
 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!-------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(ng,ng,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define (spectral) residual PV qd = (1-F)[qs-qc]:
qd=bfhi*(qs-qc)
 !Here bfhi = 1-F is a high-pass spectral filter

return
end subroutine

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (spectral) PV anomaly field is stored in qs, and the 
! residual PV needed in congen.f90 is stored in qr.

implicit none

 !Local variables:
double precision:: qa(ng,ng)

!-----------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)

 !Convert qa to spectral space as qc (note, qa is overwritten):
call ptospc(ng,ng,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define spectral PV anomaly (qs) and PV residual (qd):
qs=bflo*qs+bfhi*qc+qd
qd=qs-qc

 !Convert qd to physical space as qr (used in recontouring):
call spctop(ng,ng,qd,qr,xfactors,yfactors,xtrig,ytrig)
 !Note: qd is overwritten, but we are leaving this module next
 !      and qd will be redefined upon re-entry in subroutine init.

return
end subroutine

!=======================================================================

subroutine advance(ggen)

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method 
! (for all spectral fields, namely qs, qd, ds & gs).

! Uses an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

 !Passed variable:
logical:: ggen

 !Local variables:
integer,parameter:: niter=2

 !Spectral fields needed in time stepping:
double precision:: qsi(ng,ng),sqs(ng,ng)
double precision:: qdi(ng,ng),qdm(ng,ng),sqd(ng,ng)
double precision:: dsi(ng,ng),sds(ng,ng),nds(ng,ng)
double precision:: gsi(ng,ng),sgs(ng,ng),ngs(ng,ng)
 !Contour positions needed in time stepping:
double precision:: xqi(nptq),yqi(nptq)
 !Contour velocities:
double precision:: uq(nptq),vq(nptq)
 !Other local quantities:
double precision:: dsamp(0:nsamp2)
double precision:: xx,yy
integer:: i,iter,ix,iy
!-------------------------------------------------------------------
 !Invert PV and compute velocity at current time level, say t=t^n:
if (ggen) call inversion
 !If ggen is false, inversion was called previously at this time level.

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !Here F is a low pass filter (see spectral.f90)
qs=bflo*qs+bfhi*qc+qd
qd=bfhi*(qs-qc)

 !Compute twist parameter and save various diagnostics each time step:
call diagnose

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all contours and fields:

 !Contours:
call velint(uu,vv,uq,vq)
 !Here (uq,vq) stands for the velocity at time level n, i.e. u(x^n,t^n)
do i=1,nptq
   !Store x^n+0.5*dt*u(x^n,t^n) for efficiency in the iteration below:
  xx=xq(i)+dt2*uq(i)
  yy=yq(i)+dt2*vq(i)
  xqi(i)=oms*(xx-twopi*dble(int(xx*pinv)))
  yqi(i)=oms*(yy-twopi*dble(int(yy*pinv)))
   !Preliminary guess for x^{n+1}:
  xx=xq(i)+dt*uq(i)
  yy=yq(i)+dt*vq(i)
  xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
  yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
enddo

 !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
 !(sds,sgs) for divergence and acceleration divergence (ds,gs):
call source(sqs,sqd,sds,sgs)

 !If there are sample points for computing the frequency spectrum in
 !post-processing, save the divergence dd at these (found in source):
if (nsamp .gt. 0) then
  i=0
  do ix=1,ng,ng/nsamp
    do iy=1,ng,ng/nsamp
      i=i+1  
      dsamp(i)=dd(iy,ix)
    enddo
  enddo
  write(61,*) t,(dsamp(i),i=1,nsamp2)
endif

 !Update PV fields:
qsi=qs+dt2*sqs
qs=qs+dt*sqs
qdi=qd
qdm=qd+dt4*sqd
qd=diss*(qdm+dt4*sqd)-qdi

 !Update divergence and acceleration divergence:
dsi=ds
gsi=gs
nds=sds+dt4i*dsi
ngs=sgs+dt4i*gsi
sds=nds+sds          !2*N_tilde_delta
sgs=ngs+sgs          !2*N_tilde_gamma
ds=simp*(pdis*sds+sgs)-dsi       !simp = 1/(P*R^2-G);  pdis = P*R
gs=simp*(pdis*sgs+pgop*sds)-gsi  !pgop = P*G

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion

   !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
   !(sds,sgs) for divergence and acceleration divergence (ds,gs):
  call source(sqs,sqd,sds,sgs)

   !Interpolate gridded velocity (uu,vv) at contour nodes as (uq,vq):
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
    yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
  enddo
   !Now (xq,yq) contain a new guess for x^{n+1}.

   !Update PV fields:
  qs=qsi+dt2*sqs
  qd=diss*(qdm+dt4*sqd)-qdi

   !Update divergence and acceleration divergence:
  sds=nds+sds          !2*N_tilde_delta
  sgs=ngs+sgs          !2*N_tilde_gamma
  ds=simp*(pdis*sds+sgs)-dsi       !simp = 1/(P*R^2-G);  pdis = P*R
  gs=simp*(pdis*sgs+pgop*sds)-gsi  !pgop = P*G
enddo

 !Advance time:
t=t+dt

return
end subroutine

!=======================================================================

subroutine source(sqs,sqd,sds,sgs)

! Gets the source terms (sqs,sqd) for the PV (qs,qd) as well as those
! (sds,sgs) for divergence and acceleration divergence (ds,gs) ---
! --- all in spectral space.  

! Note that (sds,sgs) only include the nonlinear terms for a 
! semi-implicit treatment, closely analogous to that described in 
! the appendix of Mohebalhojeh & Dritschel (2004).

! The spectral fields ds, gs, qd and qs are all spectrally truncated.
! Note, hh, uu, vv, qq & zz obtained by main_invert before calling this 
! routine are all spectrally truncated.

implicit none

double precision,parameter:: toler=1.d-11
 !toler: maximum error in iteration below to find the vertically-
 !       integrated non-hydrostatic pressure \bar{P}_n (ppn below).
double precision,parameter:: w=0.75d0, wc=one-w
 !w: weight used in the pressure iteration to accelerate convergence

 !Passed variables:
double precision:: sqs(ng,ng),sqd(ng,ng),sds(ng,ng),sgs(ng,ng)

 !Local variables (physical):
double precision:: htot(ng,ng),hinv(ng,ng),hinvsq(ng,ng)
double precision:: rx(ng,ng),ry(ng,ng)
double precision:: pnx(ng,ng),pny(ng,ng),wkp(ng,ng),wkq(ng,ng)
double precision:: errpn

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng),wkf(ng,ng),wkg(ng,ng)

!---------------------------------------------------------------
 !qd source --- only NL advection term is needed:
call gradient(qd,rx,ry)
wkp=-uu*rx-vv*ry
 !Convert to spectral space:
call ptospc(ng,ng,wkp,sqd,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sqd=filt*sqd

!---------------------------------------------------------------
 !qs source --- only NL advection term is needed:
call gradient(qs,rx,ry)
wkp=-uu*rx-vv*ry
 !Convert to spectral space:
call ptospc(ng,ng,wkp,sqs,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sqs=filt*sqs

!---------------------------------------------------------------
 !Calculate the vertically-integrated non-hydrostatic pressure.

 !First form gamma_tilde = gamma + 2*(J(u,v) - delta^2):
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
 !dd contains the divergence (delta) in physical space

call jacob(uu,vv,wkp)
 !wkp contains J(u,v) in physical space

wkp=wkp-dd**2
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
wkg=gs+two*filt*wkg
call spctop(ng,ng,wkg,wkp,xfactors,yfactors,xtrig,ytrig)
 !wkp now contains gamma_tilde in physical space (de-aliased)

 !Multiply next by (1+h) and re-store in wkg (spectral) as the 
 !fixed rhs in the pressure iteration below:
htot=one+hh
wkp=htot*wkp
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
 !wkg is not de-aliased by the prop operator below takes care of this.

 !Next calculate wkq = 3/H^2*(1/(1+h)^2 - 1) needed for the pressure
 !iteration below:
hinv=one/htot
call dealias(hinv)
hinvsq=hinv**2
call dealias(hinvsq)
wkq=hbsq3i*(hinvsq-one)

 !Calculate also (1+h)^{-1}*grad(h) and store in rx & ry:
wkp=hh
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,rx,ry)
rx=hinv*rx
call dealias(rx)
ry=hinv*ry
call dealias(ry)

 !Now iterate to find \bar{P}_n (in ppn) starting from the guess
 !ppn = (grad^2 - 3/H^2)^{-1}((1+h)*gamma_tilde):
wka=prop*wkg
call spctop(ng,ng,wka,ppn,xfactors,yfactors,xtrig,ytrig)
errpn=two*toler
do while (errpn .gt. toler)
  wkp=ppn
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call gradient(wka,pnx,pny)
  wkp=rx*pnx+ry*pny+wkq*ppn
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  wka=prop*(wkg+wkb)
  call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
  errpn=sqrt(dsumi*sum((wkp-ppn)**2))
  ppn=w*wkp+wc*ppn
enddo

 !Store rhs of pressure iteration, minus gs, in wkb for use in divergence
 !tendency below:
wkb=wkg-gs+wkb

 !Compute J(ppn,hinv) and store spectral version in wkf (not de-aliased):
call jacob(ppn,hinv,wkp)
call ptospc(ng,ng,wkp,wkf,xfactors,yfactors,xtrig,ytrig)
 !wkf is needed below for the gamma source.

!---------------------------------------------------------------
 !Nonlinear part of ds source

 !Compute div(delta*u,delta*v) and store in wka (spectral):
wkp=dd*uu
wkq=dd*vv
call divs(wkp,wkq,wka)

 !Compute delta^2 and store in spectral space as wkg:
wkp=dd**2
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)

 !Store (spectral) 2*delta^2 - div(delta*u,delta*v) in wka to free up wkg:
wka=two*wkg-wka

 !Compute 1/(1+h)^3 -> wkp and de-alias:
wkp=hinv*hinvsq
call dealias(wkp)
 !Form \bar{P}_n*(1/(1+h)^3 - 1):
wkp=ppn*(wkp-one)
 !Transform to spectral space as wkg:
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)

 !Add everything up to define delta_t - pope*gamma (spectral, filtered)
sds=filt*(wka-hbsq3i*wkg+pope*wkb)
 !Here pope = (1 - (H^2/3)*grad^2)^{-1} (spectral)

!---------------------------------------------------------------
 !Nonlinear part of gs source --- f*N_zeta - c^2*Lap(N_h) where
 !N_zeta and N_h are the nonlinear parts of the vorticity and
 !depth sources:

 !Obtain N_h:
wkp=hh*uu
wkq=hh*vv
 !Compute div(h*u,h*v) = -N_h spectrally and filter for use below:
call divs(wkp,wkq,wka)
 !For use below, calculate -c^2*Lap(N_h):
wka=c2g2*wka

 !Obtain N_zeta:
wkp=zz*uu
wkq=zz*vv
 !Compute div(zeta*u,zeta*v) spectrally:
call divs(wkp,wkq,wkb)

 !Add everything up to define S_gamma:
sgs=wka+cof*filt*(wkf-wkb)
 !wkf contains J(h,2*P_0/3)/(1+h) in spectral space, part of N_zeta.

return
end subroutine

!=======================================================================

subroutine nhpsolve

! Finds the NH pressure (used only diagnostically in data writes).
! The physical space diverence is assumed to be available in dd.
  
implicit none

double precision,parameter:: toler=1.d-11
 !toler: maximum error in iteration below to find the vertically-
 !       integrated non-hydrostatic pressure \bar{P}_n (ppn below).
double precision,parameter:: w=0.75d0, wc=one-w
 !w: weight used in the pressure iteration to accelerate convergence

 !Local variables (physical):
double precision:: htot(ng,ng),hinv(ng,ng),hinvsq(ng,ng)
double precision:: rx(ng,ng),ry(ng,ng)
double precision:: pnx(ng,ng),pny(ng,ng),wkp(ng,ng),wkq(ng,ng)
double precision:: errpn

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng),wkg(ng,ng)

!---------------------------------------------------------------
 !Calculate the vertically-integrated non-hydrostatic pressure.

 !First form gamma_tilde = gamma + 2*(J(u,v) - delta^2):
call jacob(uu,vv,wkp)
 !wkp contains J(u,v) in physical space

wkp=wkp-dd**2
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
wkg=gs+two*filt*wkg
call spctop(ng,ng,wkg,wkp,xfactors,yfactors,xtrig,ytrig)
 !wkp now contains gamma_tilde in physical space (de-aliased)

 !Multiply next by (1+h) and re-store in wkg (spectral) as the 
 !fixed rhs in the pressure iteration below:
htot=one+hh
wkp=htot*wkp
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
 !wkg is not de-aliased by the prop operator below takes care of this.

 !Next calculate wkq = 3/H^2*(1/(1+h)^2 - 1) needed for the pressure
 !iteration below:
hinv=one/htot
call dealias(hinv)
hinvsq=hinv**2
call dealias(hinvsq)
wkq=hbsq3i*(hinvsq-one)

 !Calculate also (1+h)^{-1}*grad(h) and store in rx & ry:
wkp=hh
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,rx,ry)
rx=hinv*rx
call dealias(rx)
ry=hinv*ry
call dealias(ry)

 !Now iterate to find \bar{P}_n (in ppn) starting from the guess
 !ppn = (grad^2 - 3/H^2)^{-1}((1+h)*gamma_tilde):
wka=prop*wkg
call spctop(ng,ng,wka,ppn,xfactors,yfactors,xtrig,ytrig)
errpn=two*toler
do while (errpn .gt. toler)
  wkp=ppn
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call gradient(wka,pnx,pny)
  wkp=rx*pnx+ry*pny+wkq*ppn
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  wka=prop*(wkg+wkb)
  call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
  errpn=sqrt(dsumi*sum((wkp-ppn)**2))
  ppn=w*wkp+wc*ppn
enddo

return
end subroutine

!=======================================================================

subroutine inversion

! Finds the gridded dimensionless height anomaly (hh) and velocity
! field (uu,vv) from the PV contours and the (spectral) divergence (ds)
! and acceleration divergence (gs).

implicit none

 !Local variables:
double precision:: qa(ng,ng)

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is overwritten):
call ptospc(ng,ng,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Combine fields to update qt with full (spectral) field,
 !qt = F[qs-qc]+qc+qd, where F is a low pass filter:
qt=bflo*qs+bfhi*qc+qd

 !Invert PV, divergence and acceleration divergence to obtain the
 !dimensionless height anomaly and velocity field, as well as the
 !gridded PV anomaly and relative vorticity (see spectral.f90):
call main_invert(qt,ds,gs,hh,uu,vv,qq,zz)
 !Note: qt, ds & gs are in spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

return
end subroutine

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: wka(ng,ng)
double precision:: ro,fr,hmin,hmax,cfl

!----------------------------------------------------------------------
 !Compute Rossby number:
ro=maxval(abs(zz))/cof

 !Compute Froude number:
wka=uu**2+vv**2
fr=sqrt(maxval(wka/(one+hh)))/cgw

 !Compute h_min and h_max:
hmin=minval(hh)
hmax=maxval(hh)

 !Compute CFL parameter, u_max*dt/dx:
cfl=sqrt(maxval(wka))*dt*gli

 !Write data:
write(16,'(1x,f12.5,5(1x,f12.8))') t,ro,fr,hmin,hmax,cfl

 !Increment the integral of |zeta|_max:
twist=twist+dt*maxval(abs(zz))

return
end subroutine

!=======================================================================

subroutine savegrid(igrids)

! Saves PV, energy and various spectra at the desired save time

implicit none

 !Passed variable:
integer:: igrids

 !Local variables:
double precision:: htot(ng,ng),wkp(ng,ng)           !Physical
double precision:: wks(ng,ng)                       !Spectral
double precision:: spec(0:ng)                       !For 1D spectra
double precision:: ekin,epot,ediv,etot              !Energy components
integer:: k

!---------------------------------------------------------------
 !Compute energy components and total:
wks=ds
call spctop(ng,ng,wks,dd,xfactors,yfactors,xtrig,ytrig)
htot=one+hh
ekin=f12*garea*sum(htot*(uu**2+vv**2))
epot=f12*garea*csq*sum(hh**2)
wkp=htot*dd
ediv=f12*garea*hbsq3*sum(htot*wkp**2)
etot=ekin+epot+ediv

 !Write energies to ecomp.asc:
write(15,'(f9.2,5(1x,f16.9))') t,ediv,ekin,ekin+ediv,epot,etot
write(*,'(a,f9.2,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Write various gridded fields to direct access files:
write(31,rec=igrids) t,qq
write(32,rec=igrids) t,dd
wks=gs
call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
write(33,rec=igrids) t,wkp
write(34,rec=igrids) t,hh
write(35,rec=igrids) t,zz
call nhpsolve
write(36,rec=igrids) t,ppn

!---------------------------------------------------------------
 !Compute 1d spectra for various fields:
wkp=zz
call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
call spec1d(wks,spec)
spec=log10(spmf*spec+1.d-32)
write(51,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),spec(k)
enddo

call spec1d(ds,spec)
spec=log10(spmf*spec+1.d-32)
write(52,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(52,'(2(1x,f12.8))') alk(k),spec(k)
enddo

call spec1d(gs,spec)
spec=log10(spmf*spec+1.d-32)
write(53,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(53,'(2(1x,f12.8))') alk(k),spec(k)
enddo

wkp=hh
call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
call spec1d(wks,spec)
spec=log10(spmf*spec+1.d-32)
write(54,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(54,'(2(1x,f12.8))') alk(k),spec(k)
enddo

wkp=ppn
call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
call spec1d(wks,spec)
spec=log10(spmf*spec+1.d-32)
write(55,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(55,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !spmf takes into account uneven sampling of wavenumbers in each
 !shell [k-1/2,k+1/2].

 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated

 !alk(k) = log_10(k)

return
end subroutine

!=======================================================================
      
subroutine savecont(irec)

! Saves PV contours for post-processing and imaging

implicit none

 !Passed variable:
integer:: irec

 !Local variables:
double precision:: ss(ng,ng)
double precision:: qa(ng,ng)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
ss=qt-qc
call spctop(ng,ng,ss,qa,xfactors,yfactors,xtrig,ytrig)
write(83,rec=irec) t,qa

 !Save PV contours:
open(81,file='contours/qqindex'//pind,form='unformatted',status='replace')
write(81) npq(1:nq),i1q(1:nq),indq(1:nq)
close(81)

open(82,file='contours/qqnodes'//pind,form='unformatted',status='replace')
write(82) xq(1:nptq),yq(1:nptq)
close(82)

return
end subroutine

!=======================================================================

 !Main end module
end module
