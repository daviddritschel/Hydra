module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

 !Various spectral PV fields (total, contour part, residual):
double precision:: qt(ng,ng),qc(ng,ng),qd(ng,ng)

 !Spectral operators used in time stepping:
double precision:: rdis(ng,ng),simp(ng,ng),disq(ng,ng),disb(ng,ng)

 !Physical space velocity divergence:
double precision:: dd(ng,ng)

 !Logicals for saving gridded & contour data
logical:: gsave,csave

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
double precision,parameter:: qratmax=0.2d0
!      qratmax:  the maximum ratio r of the mean-square residual PV qd
!                to the mean-square PV qc in the contours
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, or when r > qratmax, 
!      the code rebuild the PV contours in a separate memory space.
double precision:: qa(ng,ng),qrat
double precision:: wka(ng,ng)
integer:: ireg

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

   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

     !Compute ratio of mean-square residual and contour PV:
    call con2grid(qa)
    wka=qd
    call spctop(ng,ng,wka,qr,xfactors,yfactors,xtrig,ytrig)
    qrat=sum(qr**2)/sum(qa**2)
    
     !Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax .or. qrat .gt. qratmax) then
       !Prepare PV residual qr for recontouring (and preserve qs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record active contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance
  
enddo
 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Save final data if not done already:
if (int(t/tgsave) .eq. igrids) call savegrid
if (int(t/tcsave) .eq. iconts) call savecont

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration

implicit none

 !Local variable:
double precision:: qa(ng,ng)

!------------------------------------------------------------------
 !Record contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!-------------------------------------------------------------
 !Logicals used for saving gridded fields and contours:
gsave=.false.
csave=.false.

 !Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(ng,ng,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define (spectral) residual PV qd = (1-F)[qs-qc]:
qd=bfhi*(qs-qc)
 !Here bfhi = 1-F is a high-pass spectral filter

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (spectral) PV anomaly field is stored in qs, and the 
! residual PV needed in congen.f90 is stored in qr.

implicit none

 !Local variables:
double precision:: qa(ng,ng)
double precision:: wka(ng,ng),qmin,ql1,ql2
integer:: ix,iy

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

!-------------------------------------------------------------------
 !Update PV contour interval for re-contouring; choose contour
 !interval based on <q^2>/<|q|> for |q| > q_rms:
wka=qs
call spctop(ng,ng,wka,qa,xfactors,yfactors,xtrig,ytrig)
qmin=sqrt(dsumi*sum(qa**2))
ql1=zero
ql2=zero
do ix=1,ng
  do iy=1,ng
    if (abs(qa(iy,ix)) .gt. qmin) then
      ql1=ql1+abs(qa(iy,ix))
      ql2=ql2+qa(iy,ix)**2
    endif
  enddo
enddo
qjump=ql2/(ql1*dble(ncont))

return
end subroutine prepare

!=======================================================================

subroutine advance

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method 
! (for all spectral fields: qs, qd, ds, gs, bxs, bys & bzs).

! Uses an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

 !Local variables:
integer,parameter:: niter=2

 !Spectral fields needed in time stepping:
double precision:: qsi(ng,ng),sqs(ng,ng)
double precision:: qdi(ng,ng),qdm(ng,ng),sqd(ng,ng)
double precision:: bxsi(ng,ng),bxsm(ng,ng),sbxs(ng,ng)
double precision:: bysi(ng,ng),bysm(ng,ng),sbys(ng,ng)
double precision:: bzsi(ng,ng),bzsm(ng,ng),sbzs(ng,ng)
double precision:: dsi(ng,ng),sds(ng,ng),nds(ng,ng)
double precision:: gsi(ng,ng),sgs(ng,ng),ngs(ng,ng)
 !Contour positions needed in time stepping:
double precision:: xqi(nptq),yqi(nptq)
 !Contour velocities:
double precision:: uq(nptq),vq(nptq)
 !Other local quantities:
double precision:: xx,yy
integer:: i,iter

!-------------------------------------------------------------------
 !Invert PV and compute velocity at current time level, say t=t^n:
call inversion
 !This also returns qc, needed immediately below.

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !Here F is a low pass filter (see spectral.f90)
qs=bflo*qs+bfhi*qc+qd
qd=bfhi*(qs-qc)

 !Adapt the time step and save various diagnostics each time step:
call adapt

 !Possibly save data (gsave & csave set by adapt):
if (gsave) call savegrid
if (csave) call savecont

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
 !(sds,sgs,sbxs,sbys,sbzs) for divergence and acceleration divergence (ds,gs):
call source(sqs,sqd,sds,sgs,sbxs,sbys,sbzs)

 !Update PV fields:
qsi=qs+dt2*sqs
qs=qs+dt*sqs
qdi=qd
qdm=qd+dt4*sqd
qd=disq*(qdm+dt4*sqd)-qdi

 !Update magnetic field:
bxsi=bxs
bxsm=bxs+dt4*sbxs
bxs=disb*(bxsm+dt4*sbxs)-bxsi
bysi=bys
bysm=bys+dt4*sbys
bys=disb*(bysm+dt4*sbys)-bysi
bzsi=bzs
bzsm=bzs+dt4*sbzs
bzs=disb*(bzsm+dt4*sbzs)-bzsi

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
   !(sds,sgs,sbxs,sbys,sbzs) for the other fields:
  call source(sqs,sqd,sds,sgs,sbxs,sbys,sbzs)

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
  qd=disq*(qdm+dt4*sqd)-qdi

   !Update magnetic field:
  bxs=disb*(bxsm+dt4*sbxs)-bxsi
  bys=disb*(bysm+dt4*sbys)-bysi
  bzs=disb*(bzsm+dt4*sbzs)-bzsi

   !Update divergence and acceleration divergence:
  sds=nds+sds          !2*N_tilde_delta
  sgs=ngs+sgs          !2*N_tilde_gamma
  ds=simp*(rdis*sds+sgs)-dsi       !simp = 1/(R^2-G);  rdis = R
  gs=simp*(rdis*sgs+opak*sds)-gsi  !opak = G
enddo

 !Advance time:
t=t+dt

return
end subroutine advance

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
call main_invert(qt,ds,gs,hh,uu,vv,zz)
 !Note: qt, ds & gs are in spectral space while 
 !      hh, uu, vv and zz are in physical space.

return
end subroutine inversion

!=======================================================================

subroutine adapt

! Adapts the time step, computes the twist parameter (the time integral
! of |zeta|_max) and various quantities every time step to monitor the
! flow evolution, and updates spectral operators needed in the time
! stepping of delta and gamma_l.

implicit none

 !Local variables:
double precision,parameter:: dtfac=pi/10.d0, cflmax=0.75d0
! The time step dt is chosen to be less than or equal to the minimum
! of dtfac/max(|zeta|_max) and cflmax*dx/max(|u|_max,|b|_max,c_gw)
! where dx is the grid spacing, u is the horizontal velocity field,
! b is the horizontal magnetic field, and c_gw is the short-scale
! inertia-gravity wave speed.

double precision:: wka(ng,ng),bx(ng,ng),by(ng,ng)
double precision:: bmax,umax,zmax,zrms,ro,fr,hmin,hmax
integer:: itime

!----------------------------------------------------------------------
 !Get horizontal magnetic field in physical space:
wka=bxs
call spctop(ng,ng,wka,bx,xfactors,yfactors,xtrig,ytrig)
wka=bys
call spctop(ng,ng,wka,by,xfactors,yfactors,xtrig,ytrig)

 !Maximum horizontal magnetic field:
bmax=sqrt(maxval(bx**2+by**2))

 !Maximum horizontal velocity:
wka=uu**2+vv**2
umax=sqrt(maxval(wka))

 !Froude number:
fr=sqrt(maxval(wka/(one+hh)))/cgw

 !Maximum vorticity:
zmax=maxval(abs(zz))

 !R.m.s. vorticity:
zrms=sqrt(dsumi*sum(zz**2))

 !Rossby number:
ro=zmax/cof

 !Min/max dimensionless height anomaly:
hmin=minval(hh)
hmax=maxval(hh)

 !Write data:
write(16,'(1x,f12.5,4(1x,f12.8))') t,ro,fr,hmin,hmax
write(17,'(1x,f12.5,4(1x,f12.8))') t,bmax,umax,zmax,zrms

!---------------------------------------------------------------------
 !Compute new time step:
dt=min(gl*cflmax/max(umax,bmax,cgw),dtfac/(zmax+small),dtmax)
 !Note, dtmax is specified in constants.f90

 !Time step parameters:
dt2=dt*f12
dt4=dt*f14
dt2i=one/dt2
dt4i=one/dt4

 !Update operators needed in time stepping:
rdis=dt2i+(cof+zrms)*qdis
simp=filt/(rdis**2-opak)
disq=two/(dt2*rdis)
disb=two/(one+dt2*bdis)
rdis=filt*rdis

 !Increment the integral of |zeta|_max:
twist=twist+dt*zmax

!---------------------------------------------------------------------
 !Set flags to save data:
gsave=(int(t/tgsave) .eq. igrids)
 !Gridded data will be saved at time t if gsave is true.
csave=(int(t/tcsave) .eq. iconts)
 !Contour data will be saved at time t if csave is true.

return
end subroutine adapt

!=======================================================================

subroutine source(sqs,sqd,sds,sgs,sbxs,sbys,sbzs)

! Gets the source terms (sqs,sqd) for the PV (qs,qd) as well as those
! (sds,sgs,sbxs,sbys,sbzs) for divergence ds, acceleration divergence gs
! and the magnetic field (bxs,bys,bzs) --- all in spectral space.
! Note that (sds,sgs) only include the nonlinear terms for a semi-implicit
! treatment, closely analogous to that described in the appendix of
! Mohebalhojeh & Dritschel (2004).

! The fields ds, gs, qd, qs, bxs, bys & bzs are all spectrally truncated.
! Note, hh, uu, vv & zz obtained by main_invert before calling this 
! routine are all spectrally truncated.

implicit none

 !Passed variables:
double precision:: sqs(ng,ng),sqd(ng,ng),sds(ng,ng),sgs(ng,ng)
double precision:: sbxs(ng,ng),sbys(ng,ng),sbzs(ng,ng)

 !Local variables (physical):
double precision:: wkp(ng,ng),wkq(ng,ng)
double precision:: bx(ng,ng),by(ng,ng),bz(ng,ng),jz(ng,ng)
double precision:: ux(ng,ng),uy(ng,ng),vx(ng,ng),vy(ng,ng)

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng)
double precision:: fbx(ng,ng),fby(ng,ng),sql(ng,ng)

!---------------------------------------------------------------
 !Prepare for source calculations:

 !Get divergence in physical space as dd:
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)

 !Get (B_x,B_y,b) in physical space as (bx,by,bz):
wka=bxs
call spctop(ng,ng,wka,bx,xfactors,yfactors,xtrig,ytrig)
wka=bys
call spctop(ng,ng,wka,by,xfactors,yfactors,xtrig,ytrig)
wka=bzs
call spctop(ng,ng,wka,bz,xfactors,yfactors,xtrig,ytrig)

 !Get j_z = dB_y/dx - dB_x/dy in physical space as jz:
call xderiv(ng,ng,hrkx,bys,wka)
call yderiv(ng,ng,hrky,bxs,wkb)
wka=wka-wkb
call spctop(ng,ng,wka,jz,xfactors,yfactors,xtrig,ytrig)

 !Get horizontal Lorentz force in spectral space as (fbx,fby):
call gradient(bzs,wkp,wkq) !wkp = db/dx & wkq = db/dy (physical)
wkp=-jz*by-bz*wkp
wkq= jz*bx-bz*wkq
call ptospc(ng,ng,wkp,fbx,xfactors,yfactors,xtrig,ytrig)
call ptospc(ng,ng,wkq,fby,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
fbx=filt*fbx
fby=filt*fby

 !Get S_q = curl(fbx,fby) - div(q_l*(u,v)) in spectral space as sql:
call xderiv(ng,ng,hrkx,fby,wka) !wka = dfby/dx (spectral)
call yderiv(ng,ng,hrky,fbx,wkb) !wkb = dfbx/dy (spectral)
sql=qt !qt contains q_l is spectral space
call spctop(ng,ng,sql,wkq,xfactors,yfactors,xtrig,ytrig)
 !wkq contains q_l in physical space
wkp=wkq*uu
wkq=wkq*vv
call divs(wkp,wkq,sql) !sql contains div(q_l*(u,v)) in spectral space
sql=wka-wkb-filt*sql

!---------------------------------------------------------------
 !qs source --- only NL advection term is needed:

call gradient(qs,wkp,wkq)
wkp=-uu*wkp-vv*wkq
 !Convert to spectral space:
call ptospc(ng,ng,wkp,sqs,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sqs=filt*sqs

!---------------------------------------------------------------
 !qd source --- S_q + (u,v)*grad{q_l - q_d}:

wkb=qt-qd
call gradient(wkb,wkp,wkq) !(wkp,wkq) is grad{q_l - q_d} (physical)
wkp=uu*wkp+vv*wkq
call ptospc(ng,ng,wkp,sqd,xfactors,yfactors,xtrig,ytrig)
 !Add S_q to complete qd source and apply de-aliasing filter:
sqd=sql+filt*sqd

!---------------------------------------------------------------
 !Nonlinear part of ds source, 2*J(u,v) + div(fbx,fby):

wkp=uu
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,ux,uy) !ux = du/dx & uy = du/dy (physical)
vx=zz+uy !vx = dv/dx using definition of vorticity
vy=dd-ux !vy = dv/dy using definition of diverence

 !Form J(u,v):
wkp=ux*vy-uy*vx
 !Convert J(u,v) to spectral space as sds:
call ptospc(ng,ng,wkp,sds,xfactors,yfactors,xtrig,ytrig)

 !Form 2*J(u,v) + div(fbx,fby) in spectral space as sds:
call xderiv(ng,ng,hrkx,fbx,wka) !wka = dfbx/dx (spectral)
call yderiv(ng,ng,hrky,fby,wkb) !wkb = dfby/dy (spectral)
sds=two*filt*sds+wka+wkb

!---------------------------------------------------------------
 !Nonlinear part of gs source --- f*S_q - G{div(h_tilde*(u,v))}
 !where S_q = dq_l/dt and G = c^2*grad^2 - f^2 is the spectral
 !gravity-wave operator:

 !Compute div(h_tilde*(u,v)) spectrally (put into wka):
wkp=hh*uu
wkq=hh*vv
call divs(wkp,wkq,wka)

 !Add everything up to define S_gamma:
sgs=gwop*wka+cof*sql
 !gwop = G (filtered)

!---------------------------------------------------------------
 !bxs source (except diffusion), (B_x,B_y)*grad(u)-(u,v)*grad(B_x):

call gradient(bxs,wkp,wkq)
 !(wkp,wkq) = grad{B_x} in physical space
wkp=bx*ux+by*uy-wkp*uu-wkq*vv
call ptospc(ng,ng,wkp,sbxs,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sbxs=filt*sbxs

!---------------------------------------------------------------
 !bys source (except diffusion), (B_x,B_y)*grad(v)-(u,v)*grad(B_y):

call gradient(bys,wkp,wkq)
 !(wkp,wkq) = grad{B_y} in physical space
wkp=bx*vx+by*vy-wkp*uu-wkq*vv
call ptospc(ng,ng,wkp,sbys,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sbys=filt*sbys

!---------------------------------------------------------------
 !bzs source (except diffusion), -div(b*u,b*v):

wkp=bz*uu
wkq=bz*vv
call divs(wkp,wkq,sbzs)
 !Apply de-aliasing filter:
sbzs=-filt*sbzs

return
end subroutine source

!=======================================================================

subroutine savegrid

! Saves PV, energy, various spectra and energy dissipation rates
! at the desired save time

implicit none

 !Local variables:
double precision:: htot(ng,ng),bx(ng,ng),by(ng,ng),bz(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng)
double precision:: spec(0:ng),tspec(0:ng)
double precision:: ekin,epot,emag,etot
double precision:: dbh,dbz,deta
integer:: k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Update fields:
call inversion

 !Compute energy components and total:
wkb=bxs
call spctop(ng,ng,wkb,bx,xfactors,yfactors,xtrig,ytrig)
wkb=bys
call spctop(ng,ng,wkb,by,xfactors,yfactors,xtrig,ytrig)
wkb=bzs
call spctop(ng,ng,wkb,bz,xfactors,yfactors,xtrig,ytrig)
htot=one+hh
ekin=f12*garea*sum(htot*(uu**2+vv**2))
epot=f12*garea*csq*sum(hh**2)
emag=f12*garea*sum(htot*(bx**2+by**2+bz**2))
etot=ekin+epot+emag

 !Write energies to ecomp.asc:
write(15,'(f13.6,4(1x,f16.9))') t,ekin,epot,emag,etot
write(*,'(a,f13.6,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Write various gridded fields to direct access files
 !(the magnetic quantities are written when computing spectra):
wkb=ds
call spctop(ng,ng,wkb,dd,xfactors,yfactors,xtrig,ytrig)
write(32,rec=igrids) t,dd
wkb=gs
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
write(33,rec=igrids) t,wkp
write(34,rec=igrids) t,hh
write(35,rec=igrids) t,zz

!---------------------------------------------------------------
 !Compute vorticity spectrum:
wkp=zz
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
call spec1d(wkb,spec)
spec=log10(spmf*spec+1.d-32)
write(51,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !Compute velocity divergence spectrum:
call spec1d(ds,spec)
spec=log10(spmf*spec+1.d-32)
write(52,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(52,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !Compute acceleration divergence spectrum:
call spec1d(gs,spec)
spec=log10(spmf*spec+1.d-32)
write(53,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(53,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !Compute dimensionless height anomaly spectrum:
wkp=hh
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
call spec1d(wkb,spec)
spec=log10(spmf*spec+1.d-32)
write(54,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(54,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !Get j_z = dB_y/dx - dB_x/dy in physical space as jz:
call xderiv(ng,ng,hrkx,bys,wka)
call yderiv(ng,ng,hrky,bxs,wkb)
wkb=wka-wkb
 !Compute its spectrum:
call spec1d(wkb,spec)
spec=log10(spmf*spec+1.d-32)
write(55,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(55,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !Convert j_z to physical space and write:
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
write(36,rec=igrids) t,wkp

 !Compute and write -div((1+h)*(B_x,B_y):
wkp=(one+hh)*bx
wkq=(one+hh)*by
call divs(wkp,wkq,wka)
wka=filt*wka
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
wkq=wkp !Stores div((1+h)*(B_x,B_y) for computing mag. diss. below
wkp=-wkp
write(37,rec=igrids) t,wkp

 !Compute lower surface B_z (b) spectrum:
call spec1d(bzs,spec)
spec=log10(spmf*spec+1.d-32)
write(56,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(56,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !Write gridded b:
write(38,rec=igrids) t,bz

 !spmf takes into account uneven sampling of wavenumbers in each
 !shell [k-1/2,k+1/2].

 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated

 !alk(k) = log_10(k)

!-----------------------------------------------------------
 !Compute various contributions to the magnetic dissipation:
wkp=uu*bx+vv*by
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
dbh=garea*sum(wkq*wkp)
 !Above, wkq = div((1+h)*(B_x,B_y); below it is redefined:
wkq=one+hh
wkp=bz**2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
dbz=garea*sum(wkq*wkp*dd)
 !Below, bdis = eta*(k^2+l^2) in spectal space
wka=bdis*bx
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
wkp=wkq*bx*wkp
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
deta=sum(wkp)
wka=bdis*by
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
wkp=wkq*by*wkp
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
deta=deta+sum(wkp)
wka=bdis*bz
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
wkp=wkq*bz*wkp
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
deta=deta+sum(wkp)
deta=garea*deta
write(18,'(1x,f12.5,4(1x,f14.11))') t,dbh,dbz,deta,dbh+dbz+deta

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont

! Saves PV contours for post-processing and imaging

implicit none

 !Local variables:
double precision:: ss(ng,ng)
double precision:: qa(ng,ng)
character(len=3):: pind

!---------------------------------------------------------------
 !Increment counter for direct file access:
iconts=iconts+1

write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') iconts-1

 !Write contours to the contours subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
ss=qt-qc
call spctop(ng,ng,ss,qa,xfactors,yfactors,xtrig,ytrig)
write(83,rec=iconts) t,qa

 !Save PV contours:
open(81,file='contours/qqindex'//pind,form='unformatted',status='replace')
write(81) npq(1:nq),i1q(1:nq),indq(1:nq)
close(81)

open(82,file='contours/qqnodes'//pind,form='unformatted',status='replace')
write(82) xq(1:nptq),yq(1:nptq)
close(82)

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
