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
 !Physical fields:
double precision:: qx(ng,ng),qy(ng,ng)
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
ds=simp*(rdis*sds+sgs)-dsi       !simp = 1/(R^2-G);  rdis = R
gs=simp*(rdis*sgs+opak*sds)-gsi  !opak = G

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
  ds=simp*(rdis*sds+sgs)-dsi       !simp = 1/(R^2-G);  rdis = R
  gs=simp*(rdis*sgs+opak*sds)-gsi  !pgop = G
enddo

 !Advance time:
t=t+dt

return
end subroutine

!=======================================================================

subroutine source(sqs,sqd,sds,sgs)

! Gets the source terms (sqs,sqd) for the PV (qs,qd) as well as those
! (sds,sgs) for divergence and acceleration divergence (ds,gs) ---
! --- all in spectral space.  Note that (sds,sgs) only include the
! nonlinear terms for a semi-implicit treatment, closely analogous
! to that described in the appendix of Mohebalhojeh & Dritschel (2004).

! The spectral fields ds, gs, qd and qs are all spectrally truncated.
! Note, hh, uu, vv & zz obtained by main_invert before calling this 
! routine are all spectrally truncated.

implicit none

 !Passed variables:
double precision:: sqs(ng,ng),sqd(ng,ng),sds(ng,ng),sgs(ng,ng)

 !Local variables (physical):
double precision:: wkp(ng,ng),wkq(ng,ng)

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng)

!---------------------------------------------------------------
 !qd source --- only NL advection term is needed:
call gradient(qd,wkp,wkq)
wkp=-uu*wkp-vv*wkq
 !Convert to spectral space:
call ptospc(ng,ng,wkp,sqd,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sqd=filt*sqd

!---------------------------------------------------------------
 !qs source --- only NL advection term is needed:
call gradient(qs,wkp,wkq)
wkp=-uu*wkp-vv*wkq
 !Convert to spectral space:
call ptospc(ng,ng,wkp,sqs,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sqs=filt*sqs

!---------------------------------------------------------------
 !Nonlinear part of ds source --- 2J(u,v) - div(delta*(u,v)):
call jacob(uu,vv,wkp)
 !Convert J(u,v) (in wkp) to spectral space as wkb:
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)

 !Compute div(delta*u,delta*v) (to be put into wka, in spectral space):
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
 !dd contains the divergence in physical space
wkp=dd*uu
wkq=dd*vv
 !Compute spectral divergence from physical fields:
call divs(wkp,wkq,wka)

 !Add everything up to define delta_t - gamma (spectral, filtered)
sds=filt*(two*wkb-wka)

!---------------------------------------------------------------
 !Nonlinear part of gs source --- fN_zeta - c^2*Lap(N_h) where
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
sgs=wka-cof*filt*wkb

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
double precision:: wkp(ng,ng)          !Physical
double precision:: wks(ng,ng)          !Spectral
double precision:: spec(0:ng)          !For 1D spectra
double precision:: ekin,epot,etot      !Energy components
integer:: k

!---------------------------------------------------------------
 !Compute energy components and total:
ekin=f12*garea*sum((one+hh)*(uu**2+vv**2))
epot=f12*garea*csq*sum(hh**2)
etot=ekin+epot

 !Write energies to ecomp.asc:
write(15,'(f13.6,5(1x,f16.9))') t,zero,ekin,ekin,epot,etot
write(*,'(a,f13.6,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Write various gridded fields to direct access files:
write(31,rec=igrids) t,qq
wks=ds
call spctop(ng,ng,wks,dd,xfactors,yfactors,xtrig,ytrig)
write(32,rec=igrids) t,dd
wks=gs
call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
write(33,rec=igrids) t,wkp
write(34,rec=igrids) t,hh
write(35,rec=igrids) t,zz

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

 !Write contours to the contours subdirectory:
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
