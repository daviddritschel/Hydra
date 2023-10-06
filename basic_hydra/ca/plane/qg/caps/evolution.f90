module evolution

! Module contains subroutines to evolve PV contours and all fields 
! according to the algorithm detailed in caps.f90.

use common

implicit none

 !Velocity field:
double precision:: uu(ny,nx),vv(ny,nx)

 !Spectral fields (note array order):
double precision:: pp(nx,ny),qq(nx,ny),qc(nx,ny),qd(nx,ny)
double precision:: qspre(nx,ny)
double precision:: emq(nx,ny),epq(nx,ny)

 !Energies:
double precision:: ekinpre,epotpre

 !Logicals to indicate presence of contours or data saves:
logical:: contex,gsave,csave

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
  call advance
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Save final data if not already done:
gsave=abs(tinit+tgsave*dble(igrec)-tsim) .lt. 1.d-6*dt
csave=abs(tinit+tcsave*dble(icrec)-tsim) .lt. 1.d-6*dt
if (gsave .or. csave) call inversion
if (gsave) call savegrid
if (csave) call savecont

return
end subroutine

!=======================================================================

subroutine init

! Initialises quantities needed for normal time integration following 
! contour regeneration

implicit none

 !Local variable:
double precision:: qa(ny,nx)

!------------------------------------------------------------------
 !Logical to indicate presence of contours:
contex=(nptq .gt. 0)

 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i9,1x,i10)') t,nq,nptq

!-------------------------------------------------------------
 !Logicals used for saving gridded fields and contours:
gsave=.false.
csave=.false.

 !Convert PV contours to gridded values (qa):
call con2grid(qa)

 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define residual PV qd = (1-F)[qs-qc], and copy current gridded 
 !field (qs) for use in time interpolation:
qd=fhi*(qs-qc)
qspre=qs 
 !Here fhi = 1-F is a high-pass spectral filter

return
end subroutine

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current PV anomaly field is stored in qs, the residual PV needed 
! in congen.f90 is stored in qr, and (if present) tracer contours are
! stored in tracer.bin

implicit none

 !Local variable:
double precision:: qa(ny,nx)

!-----------------------------------------------------------------
 !Convert PV contours to gridded values (qa):
call con2grid(qa)

 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Put current PV anomaly (spectral in qs) and define residual (qd):
qs=flo*(qs-qc)+qc+qd
qd=qs-qc

 !Convert qd to physical space as qr (used in recontouring):
call spctop(nx,ny,qd,qr,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!=======================================================================

subroutine advance

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method (for all
! gridded fields, i.e. qs & qd).

! *** Uses a 4th-order Runge-Kutta method ***

implicit none

 !Local variables:

 !Spectral fields needed in Runge-Kutta time stepping (note array order):
double precision:: qsi(nx,ny),qsf(nx,ny),sqs(nx,ny)
double precision:: qdi(nx,ny),qdf(nx,ny),sqd(nx,ny)
 !Contour positions needed in Runge-Kutta time stepping:
double precision:: xqi(nptq),yqi(nptq),xqf(nptq),yqf(nptq)
 !Contour velocities:
double precision:: uq(nptq),vq(nptq)
 !Other local quantities:
double precision:: xx,yy
integer:: i

!-------------------------------------------------------------------
 !RK4 predictor step to time t0 + dt/2:

 !Invert PV and compute velocity:
call inversion

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !after qs is reset; here F is a low pass filter (see spectral.f90)
qs=flo*(qs-qc)+qc+qd
qd=fhi*(qs-qc)

 !Possibly save data (gsave & csave set by adapt in the previous time step):
if (gsave) call savegrid
if (csave) call savecont

 !Adjust timestep (dt) on maximum vorticity magnitude or CFL:
call adapt

 !Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,0)

if (contex) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xqi(i)=xq(i)
    yqi(i)=yq(i)
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    xx=xqi(i)+dt6*uq(i)
    yy=yqi(i)+dt6*vq(i)
    xqf(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqf(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

qsi=qs
qs=qsi+dt2*sqs
qsf=qsi+dt6*sqs
qdi=qd
qd=emq*(qdi+dt2*sqd)
qdf=qdi+dt6*sqd

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt/2:
t=t+dt2

 !Invert PV and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,1)

if (contex) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    xx=xqf(i)+dt3*uq(i)
    yy=yqf(i)+dt3*vq(i)
    xqf(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqf(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

qs=qsi+dt2*sqs
qsf=qsf+dt3*sqs
qd=emq*(qdi+dt2*sqd)
qdf=qdf+dt3*sqd

!------------------------------------------------------------------
 !RK4 predictor step at time t0 + dt:

 !Invert PV and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,1)

if (contex) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt*uq(i)
    yy=yqi(i)+dt*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    xx=xqf(i)+dt3*uq(i)
    yy=yqf(i)+dt3*vq(i)
    xqf(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqf(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

qs=qsi+dt*sqs
qsf=qsf+dt3*sqs
emq=emq**2
qd=emq*(qdi+dt*sqd)
qdf=qdf+dt3*sqd

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt:
t=t+dt2

 !Invert PV and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,2)

if (contex) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqf(i)+dt6*uq(i)
    yy=yqf(i)+dt6*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

qs=qsf+dt6*sqs
qd=emq*(qdf+dt6*sqd)

 !Add any stochastic forcing (if present) to qd:
if (stoch) call pvforce

 !Add any narrow-band spectral forcing to qd:
if (forcing) then
  call nb_forcing(dt)
  qd=qd+dqdt
endif

return
end subroutine

!=======================================================================

subroutine pvforce

! Adds stochastic forcing to PV in the array qd at the end of a time step

implicit none

 !Local variables:
double precision:: wkp(ny,nx) !Physical work array
double precision:: wks(nx,ny) !Spectral work array
double precision:: dnx,dny,svor,xx,pxc,px,yy,pyc,py
double precision:: theta,cth,sth
integer:: nvor,k,ix,iy,ix0,ix1,iy0,iy1

!-----------------------------------------------------------------
 !Initialise array to uptake added PV:
do ix=1,nx
  do iy=1,ny
    wkp(iy,ix)=zero
  enddo
enddo

 !Add point vortices here to force the flow (they are converted to 
 !gridded PV values and added to qd):
if (ivor .eq. 1) then
   !Add nvor vortices (as +/- arbitrarily separated pairs):
  nvor=2*nint(f12*(dnvor*t-totnvor))

   !The vortices are place randomly in the domain:
  dnx=dble(nx)
  dny=dble(ny)
  svor=-vorvor
  do k=1,nvor
    svor=-svor

     !Divide each vortex's circulation among the corners
     !of the local grid box (inverse bi-linear interpolation);
    xx=dnx*rand(0)
    ix0=1+int(xx)
    pxc=dble(ix0)-xx
    px=one-pxc
    ix1=ixp(ix0)

    yy=dny*rand(0)
    iy0=1+int(yy)
    pyc=dble(iy0)-yy
    py=one-pyc
    iy1=iyp(iy0)

    wkp(iy0,ix0)=wkp(iy0,ix0)+svor*pyc*pxc
    wkp(iy0,ix1)=wkp(iy0,ix1)+svor*pyc*px
    wkp(iy1,ix0)=wkp(iy1,ix0)+svor*py*pxc
    wkp(iy1,ix1)=wkp(iy1,ix1)+svor*py*px
  enddo

else
   !Add nvor dipoles:
  nvor=nint(dnvor*t-totnvor)

   !The dipoles are place randomly in the domain with random
   !orientation:
  dnx=dble(nx)
  dny=dble(ny)
  do k=1,nvor

    theta=twopi*rand(0)
    cth=cos(theta)
    sth=sin(theta)

    xx=dnx*rand(0)
    ix0=1+int(xx)
    pxc=dble(ix0)-xx
    px=one-pxc
    ix1=ixp(ix0)

    yy=dny*rand(0)
    iy0=1+int(yy)
    pyc=dble(iy0)-yy
    py=one-pyc
    iy1=iyp(iy0)

    wkp(iy0,ix0)=wkp(iy0,ix0)-vorvor*(pyc*cth+pxc*sth)
    wkp(iy0,ix1)=wkp(iy0,ix1)+vorvor*(pyc*cth -px*sth)
    wkp(iy1,ix0)=wkp(iy1,ix0)+vorvor*(pxc*sth -py*cth)
    wkp(iy1,ix1)=wkp(iy1,ix1)+vorvor*( px*sth +py*cth)
  enddo

endif

 !Total number of vortices/dipoles added so far:
totnvor=totnvor+dble(nvor)

 !Convert wkp to spectral space as wks:
call ptospc(nx,ny,wkp,wks,xfactors,yfactors,xtrig,ytrig)

 !Add wks to qd in spectral space:
qd=qd+wks

return
end subroutine

!=======================================================================

subroutine source(sqs,sqd,lev)

! Gets the source terms (sqs,sqd) for the PV (qs,qd)
! evolution equation (all in spectral space):

implicit none

 !Passed variables:
double precision:: sqs(nx,ny),sqd(nx,ny)
integer:: lev

 !Local variables:
double precision:: qqx(ny,nx),qqy(ny,nx)

!---------------------------------------------------------------
 !qd source:
call gradient(qd,qqx,qqy)
qqx=-uu*qqx-vv*qqy

 !Convert to spectral space:
call ptospc(nx,ny,qqx,sqd,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !qs source - only NL term is needed:
call gradient(qs,qqx,qqy)
if (beffect) then
  qqx=-uu*qqx-vv*(qqy+beta)
else
  qqx=-uu*qqx-vv*qqy
endif
 !Convert to spectral space:
call ptospc(nx,ny,qqx,sqs,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !Implement Ekman and/or thermal damping (add to sqd if present):
if (heating) then
  if (friction) then 
   !Use thermal and Ekman damping:
    sqd=sqd+therm*(pp-ppeq)-rekman*(qq+kdsq*pp)
  else
    sqd=sqd+therm*(pp-ppeq)
  endif
else
   !Only use Ekman damping
  sqd=sqd-rekman*(qq+kdsq*pp)
endif

!----------------------------------------------------------------
if (lev .eq. 0) return

 !Apply exponential integrating factors:
if (lev .eq. 1) then
  sqd=epq*sqd
else
  sqd=epq**2*sqd
endif

return
end subroutine

!=======================================================================

subroutine inversion

! Inverts Laplace's operator on PV anomaly (PV - beta*y) to obtain 
! the streamfunction (pp) ***in spectral space*** and the velocity 
! (uu,vv) = (-dpp/dy,dpp/dx) ***in physical space***

implicit none

 !Local variable:
double precision:: qa(ny,nx)

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qa)

 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Combine fields to update qq with full field,
 !qq = F[qs-qc]+qc+qd, where F is a low pass filter:
qq=flo*(qs-qc)+qc+qd

 !Invert PV to obtain velocity field: 
call main_invert(qq,uu,vv,pp)

return
end subroutine

!=======================================================================

subroutine adapt

! Adapts the time step dt to ensure that it is less than or equal to
! the minimum of dtfac/|zeta|_max and C*dx/|u|_max
! where dx is the grid spacing, u is the vector velocity field.
! C = cfl_max is specified below.

implicit none

 !Local variables:
double precision,parameter:: dtfac=pi/10.d0, cflmax=0.7d0
double precision:: zz(ny,nx) !Physical
double precision:: ss(nx,ny) !Spectral
double precision:: uumax,zzrms,zzmax,dtacc,tcont,cfl,dfac

!----------------------------------------------------------------------
 !Compute the gridded relative vorticity (zz):
if (btropic) then
   !Here kd = 0:
  ss=qq
else
   !Here kd > 0:
  ss=qq+kdsq*pp
endif
 !Above, qq = PV anomaly (q - beta*y) in spectral space.
call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)

 !Compute accurate advection time step:
zz=zz**2
zzrms=sqrt(dsumi*sum(zz))
zzmax=sqrt(maxval(zz))+small
zz=uu**2+vv**2
uumax=sqrt(maxval(zz))+small

dtacc=min(glx*cflmax/uumax,dtfac/max(zzmax,srwfm))
 !The restriction on the maximum Rossby wave frequency (srwfm)
 !ensures that the fastest Rossby wave frequency is resolved.

!---------------------------------------------------------------------
 !Choose a new time step, dt:
if (dt .gt. zero) then
  dt=min(dtacc,dtmax)
  if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
else
   !Limit max timestep to a data save time/5:
  dtmax=0.2d0*min(tgsave,tcsave)
  dt=min(dtacc,dtmax)
endif
 !Fractional time steps used in 4th-order Runge-Kutta time stepping:
dt2=dt*f12
dt3=dt*f13
dt6=dt*f16

 !Increment the integral of max|zz|:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
 !Record various diagnostics to monitor.asc:
cfl=uumax*dt/glx
write(17,'(1x,f12.5,1x,f6.4,4(1x,f12.6),1x,f5.3)') & 
     & t,cfl,f12*zzrms**2,zzrms,zzmax,uumax,twist

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
tgrid=tinit+tgsave*dble(igrec)
gsave=t+dt .ge. tgrid
if (gsave) then
   !Copy current gridded fields for use in time interpolation:
  qspre=qs
   !Compute energy:
  call energy(ekinpre,epotpre)
endif

 !Set flag to save contour data every tcsave time units:
tcont=tinit+tcsave*dble(icrec)
csave=t+dt .ge. tcont

!---------------------------------------------------------------------
 !Define spectral integrating factors used in Runge-Kutta integration:
dfac=dt2*zzrms
epq=exp(dfac*qdiss)
emq=one/epq

return
end subroutine

!=======================================================================

subroutine energy(ekin,epot)

! This routine computes the kinetic, potential and total energy from uu, 
! vv & pp (if kd > 0).

implicit none

 !Passed variables:
double precision:: ekin,epot

 !Local variables:
double precision:: ss(nx,ny) !Spectral
double precision:: zz(ny,nx) !Physical

!----------------------------------------------------------------------
 !Compute kinetic energy:
zz=uu**2+vv**2
ekin=f12*garea*sum(zz)

 !Compute potential energy:
if (btropic) then
   !kd = 0 and hence there is no potential energy:
  epot=zero
else
   !Get streamfunction pp in physical space as zz:
  ss=pp
  call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)
  zz=zz**2
  epot=f12*garea*kdsq*sum(zz)
endif

return
end subroutine

!=======================================================================

subroutine savegrid

! Saves PV, energy and various spectra at the desired save time

implicit none

 !Local variables:
double precision:: wka(ny,nx),wkb(ny,nx) !Physical
double precision:: qqs(nx,ny) !Spectral
double precision:: qspec(0:max(nx,ny))
double precision:: qql2,pt,ptc,ekin,epot
double precision:: ekinpost,epotpost,sumqspec
real:: qqr4(ny,nx),tr4
integer:: ix,iy,k

!---------------------------------------------------------------
 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Interpolate PV anomaly at save time:
qqs=pt*qspre+ptc*qq

!---------------------------------------------------------------
 !Compute kinetic and potential energies:
call energy(ekinpost,epotpost)

 !Compute time interpolated energies:
ekin=pt*ekinpre+ptc*ekinpost
epot=pt*epotpre+ptc*epotpost

 !Write energies to ene.asc:
write(15,'(f9.2,3(1x,f17.11))') tgrid,ekin,epot,ekin+epot

!---------------------------------------------------------------
 !Compute 1d PV spectrum:
call spec1d(qqs,qspec,0)
sumqspec=zero
do k=1,kmax
  sumqspec=sumqspec+qspec(k)
   !Normalise to take into account uneven sampling of wavenumbers 
   !in each shell [k-1/2,k+1/2]:
  qspec(k)=spmf(k)*qspec(k)
enddo
sumqspec=8.d0*sumqspec*dsumi

!---------------------------------------------------------------
!Write full PV to qq.r4:
call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)

tr4=real(tgrid)
do ix=1,nx
  do iy=1,ny
    qqr4(iy,ix)=real(wka(iy,ix)+bety(iy))
  enddo
enddo
write(31,rec=igrec+1) tr4,qqr4

 !Compute domain integral of (q-beta*y)^2 and write to norms.asc:
call l2norm(wka,qql2)
write(16,'(f9.2,1x,f16.9)') tgrid,qql2

write(*,'(a,f9.2,a,f13.6,a,f14.9)') ' t = ',tgrid,' <q^2> = ',qql2, &
                                    ' E_tot = ',ekin+epot

 !Write out spectrum to file:
write(51,'(f9.2,2(1x,f16.9),1x,i5)') tgrid,sumqspec,qql2,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),log10(qspec(k))
enddo
 !Note: alk(k) = log_10(k)

!--------------------------------------------------------------
 !Increment counter for direct file access:
igrec=igrec+1

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves PV contours for post-processing and imaging

implicit none

 !Local variables:
double precision:: ss(nx,ny) !Spectral
double precision:: qa(ny,nx) !Physical
real:: qdr4(ny,nx),tr4
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t

write(pind(1:3),'(i3.3)') icrec

tr4=real(t)
 !Write contours to the cont subdirectory:
write(80,'(i9,1x,i10,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
ss=qq-qc
call spctop(nx,ny,ss,qa,xfactors,yfactors,xtrig,ytrig)
qdr4=real(qa)
write(83,rec=icrec) tr4,qdr4

 !Save PV contours if any exist:
if (contex) then
  open(81,file='cont/qqindex'//pind,form='unformatted',status='replace')
  write(81) npq(1:nq),i1q(1:nq),indq(1:nq)
  close(81)

  open(82,file='cont/qqnodes'//pind,form='unformatted',status='replace')
  write(82) xq(1:nptq),yq(1:nptq)
  close(82)
endif

!--------------------------------------------------------------
 !Increment counter for naming next direct-access output files:
icrec=icrec+1

return
end subroutine

!=======================================================================

 !Main end module
end module
