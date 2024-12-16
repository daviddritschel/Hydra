module evolution

! Module contains subroutines to evolve buoyancy contours and all fields 
! according to the algorithm detailed in caps.f90.

use common

implicit none

 !Velocity field:
double precision:: uu(ny,nx),vv(ny,nx)

 !Spectral fields (note array order):
double precision:: pp(nx,ny),qq(nx,ny),qc(nx,ny),qd(nx,ny)
double precision:: emq(nx,ny),epq(nx,ny)
 !Tracer auxiliary arrays (if present):
double precision,allocatable,dimension(:,:):: emc,epc

 !Energy & enstrophy for time interpolation:
double precision:: enepre,enspre

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
do while (t .le. tfin)

   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

     !Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax) then
       !Prepare buoyancy residual qr for recontouring (and preserve qs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the buoyancy contours (surgery + node redistribution):
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

!Write final data; first invert buoyancy to compute velocity:
call inversion
call savegrid
call savecont

return
end subroutine advect

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
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!-------------------------------------------------------------
 !Logicals used for saving gridded fields and contours:
gsave=.false.
csave=.false.

 !Convert buoyancy contours to gridded values (qa):
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define residual buoyancy qd = (1-F)[qs-qc], and copy current gridded 
 !field (qs) for use in time interpolation:
qd=fhi*(qs-qc)
qs=filt*qs 
qspre=qs
 !Here fhi = 1-F is a high-pass spectral filter and filt is a 
 !de-aliasing filter (defined in spectral.f90)

if (tracer) then
   !Allocate memory for tracer (anomaly) evolution:
  allocate(emc(nx,ny),epc(nx,ny))
  cspre=cs
endif

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current buoyancy field is stored in qs, the residual buoyancy needed 
! in congen.f90 is stored in qr.

implicit none

 !Local variable:
double precision:: qa(ny,nx)

!-----------------------------------------------------------------
 !Convert buoyancy contours to gridded values (qa):
call con2grid(qa)

 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Put current buoyancy (spectral in qs) and define residual (qd):
qs=flo*(qs-qc)+qc+qd
qd=qs-qc

 !Convert qd to physical space as qr (used in recontouring):
call spctop(nx,ny,qd,qr,xfactors,yfactors,xtrig,ytrig)

if (tracer) then
   !De-allocate memory for tracer (anomaly) evolution:
  deallocate(emc,epc)
endif

return
end subroutine prepare

!=======================================================================

subroutine advance

! Advances buoyancy from time t to t+dt by a combination of contour 
! advection (for buoyancy contours) and the pseudo-spectral method (for all
! gridded fields, i.e. qs & qd).

! *** Uses a 4th-order Runge-Kutta method ***

implicit none

 !Local variables:

 !Spectral fields needed in Runge-Kutta time stepping (note array order):
double precision:: qsi(nx,ny),qsf(nx,ny),sqs(nx,ny)
double precision:: qdi(nx,ny),qdf(nx,ny),sqd(nx,ny)
double precision,allocatable,dimension(:,:):: csi,csf,scs
 !Contour positions needed in Runge-Kutta time stepping:
double precision:: xqi(nptq),yqi(nptq),xqf(nptq),yqf(nptq)
 !Contour velocities:
double precision:: uq(nptq),vq(nptq)
 !Other local quantities:
double precision:: xx,yy
integer:: i

!-------------------------------------------------------------------
 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !after qs is reset; here F is a low pass filter (see spectral.f90)
qs=flo*qs+fhi*qc+qd
qd=fhi*(qs-qc)

!------------------------------------------------------------------
 !RK4 predictor step to time t0 + dt/2:

 !Invert buoyancy and compute velocity:
call inversion

 !Possibly save data (gsave & csave set by adapt in the previous time step):
if (gsave) call savegrid
if (csave) call savecont

 !Adjust timestep (dt) on maximum vorticity magnitude or CFL:
call adapt

 !Calculate the source terms (sqs,sqd) for buoyancy (qs,qd):
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

if (tracer) then
   !Allocate memory for tracer (anomaly) evolution:
  allocate(csi(nx,ny),csf(nx,ny),scs(nx,ny))
   !Compute tracer source:
  call tracer_source(scs,0)
   !Evolve tracer:
  csi=cs
  cs=emc*(csi+dt2*scs)
  csf=csi+dt6*scs
endif

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt/2:
t=t+dt2

 !Invert buoyancy and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for buoyancy (qs,qd):
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

if (tracer) then
   !Compute tracer source:
  call tracer_source(scs,1)
   !Evolve tracer:
  cs=emc*(csi+dt2*scs)
  csf=csf+dt3*scs
endif

!------------------------------------------------------------------
 !RK4 predictor step at time t0 + dt:

 !Invert buoyancy and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for buoyancy (qs,qd):
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

if (tracer) then
   !Compute tracer source:
  call tracer_source(scs,1)
   !Evolve tracer:
  cs=emc*(csi+dt*scs)
  csf=csf+dt3*scs
endif

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt:
t=t+dt2

 !Invert buoyancy and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for buoyancy (qs,qd):
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

if (tracer) then
   !Compute tracer source:
  call tracer_source(scs,2)
   !Evolve tracer:
  cs=emc*(csf+dt6*scs)
   !De-allocate memory:
  deallocate(csi,csf,scs)
endif

return
end subroutine advance

!=======================================================================

subroutine source(sqs,sqd,lev)

! Gets the source terms (sqs,sqd) for the buoyancy (qs,qd)
! evolution equation (all in spectral space):

implicit none

 !Passed variables:
double precision:: sqs(nx,ny),sqd(nx,ny)
integer:: lev

 !Local variables:
double precision:: qqx(ny,nx),qqy(ny,nx)
double precision:: wkp(ny,nx)

!---------------------------------------------------------------
 !qd source:
call gradient(qd,qqx,qqy)
wkp=-uu*qqx-vv*qqy
 !Convert to spectral space:
call ptospc(nx,ny,wkp,sqd,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !qs source - only NL term is needed:
call gradient(qs,qqx,qqy)
wkp=-uu*qqx-vv*qqy
 !Convert to spectral space:
call ptospc(nx,ny,wkp,sqs,xfactors,yfactors,xtrig,ytrig)

!----------------------------------------------------------------
if (lev .eq. 0) then
 !Spectrally truncate sources:
  sqs=filt*sqs
  sqd=filt*sqd
 !Apply exponential integrating factors (and spectrally truncate sources):
else if (lev .eq. 1) then
  sqs=filt*sqs
  sqd= epq*sqd
else
  sqs=filt*sqs
  sqd= epq**2*sqd
endif

return
end subroutine source

!=======================================================================

subroutine tracer_source(scs,lev)

! Gets the source term (scs) for the tracer (anomaly) field (cs)
! evolution equation (all in spectral space):

implicit none

 !Passed variables:
double precision:: scs(nx,ny)
integer:: lev

 !Local variables:
double precision:: cx(ny,nx),cy(ny,nx)
double precision:: wkp(ny,nx)

!---------------------------------------------------------------
 !cs source:
call gradient(cs,cx,cy)
wkp=-uu*(dcdx+cx)-vv*(dcdy+cy)
 !Convert to spectral space:
call ptospc(nx,ny,wkp,scs,xfactors,yfactors,xtrig,ytrig)

!----------------------------------------------------------------
if (lev .eq. 0) then
 !Spectrally truncate source:
  scs=filt*scs
 !Apply exponential integrating factors (and spectrally truncate sources):
else if (lev .eq. 1) then
  scs=epc*scs
else
  scs=epc**2*scs
endif

return
end subroutine tracer_source

!=======================================================================

subroutine inversion

! Inverts SQG operator on buoyancy to obtain the streamfunction (pp)
! ***in spectral space*** and the velocity (uu,vv) = (-dpp/dy,dpp/dx) 
! ***in physical space***

implicit none

 !Local variable:
double precision:: qa(ny,nx)

!------------------------------------------------------------
 !Call con2grid to get updated contour buoyancy (qc):
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Combine fields to update qq with full field,
 !qq = F[qs-qc]+qc+qd, where F is a low pass filter:
qq=flo*qs+fhi*qc+qd

 !Invert buoyancy to obtain velocity field: 
call main_invert(qq,uu,vv,pp)

return
end subroutine inversion

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
double precision:: umax,zzrms,zzmax,dtacc,tcont,cfl,dfac,eif
integer:: itime

!----------------------------------------------------------------------
 !Compute the gridded relative vorticity (zz):
ss=qq*vorop
 !vorop = |k|, and is defined in spectral.f90
 !Above, qq = buoyancy in spectral space.
call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)

 !Compute accurate advection time step:
umax=sqrt(maxval(uu**2+vv**2))
zzmax=maxval(abs(zz))
zzrms=sqrt(dsumi*sum(zz**2))
dtacc=min(glx*cflmax/umax,dtfac/zzmax)

!---------------------------------------------------------------------
 !Choose a new time step, dt:
if (dt .gt. zero) then
  dt=min(dtacc,dtmax)
  if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
else
   !Limit max timestep to a data save time
  dtmax=min(tgsave,tcsave)
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
cfl=umax*dt/glx
write(17,'(1x,f12.5,1x,f12.7,1x,f14.7,1x,f12.7)') t,zzrms,zzmax,umax

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
tgrid=tgsave*dble(itime)+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  gsave=.true.

   !Copy current gridded fields for use in time interpolation:
  qspre=qs

   !Do same for the tracer (anomaly), if present:
  if (tracer) cspre=cs

   !Compute energy & enstrophy:
  call energy(enepre,enspre)
endif

 !Set flag to save contour data every tcsave time units:
itime=int((t+dt)/tcsave)
tcont=tcsave*dble(itime)
if (t .le. tcont .and. t+dt .gt. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  csave=(sign(1,2*itime-1) .gt. 0)
   !This construction avoids saving the contours at t = 0.
endif

!---------------------------------------------------------------------
 !Define spectral integrating factors used in Runge-Kutta integration:
dfac=dt2*zzrms
ss=exp(dfac*qdiss)
epq=ss*filt
emq=one/ss

 !Possibly include factors for a tracer as well:
if (tracer) then
  ss=exp(dt2*tdiss)
  epc=ss*filt
  emc=one/ss
endif

return
end subroutine adapt

!=======================================================================

subroutine energy(ene,ens)

! This routine computes the total energy & enstrophy

implicit none

 !Passed variables:
double precision:: ene,ens

 !Local variables:
double precision:: wka(ny,nx) !Physical
double precision:: qqs(nx,ny) !Spectral

!----------------------------------------------------------------------
 !Compute energy & enstrophy:
ene=f12*garea*sum(uu**2+vv**2)

qqs=qq
call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)
ens=f12*garea*sum(wka**2)

return
end subroutine energy

!=======================================================================

subroutine savegrid

! Saves buoyancy, energy and various spectra at the desired save time

implicit none

 !Local variables:
double precision:: wka(ny,nx),wkb(ny,nx) !Physical
double precision:: qqs(nx,ny),zzs(nx,ny) !Spectral
double precision:: spec(0:max(nx,ny))
double precision:: pt,ptc,ene,ens
double precision:: enepost,enspost
integer:: k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Interpolate buoyancy at save time:
qqs=pt*qspre+ptc*qq

 !Compute vertical vorticity:
zzs=vorop*qqs

!---------------------------------------------------------------
 !Compute energy:
call energy(enepost,enspost)

 !Compute time interpolated energy & enstrophy:
ene=pt*enepre+ptc*enepost
ens=pt*enspre+ptc*enspost

 !Write energy & enstrophy to ene-ens.asc:
write(15,'(f7.2,2(1x,f16.9))') tgrid,ene,ens

!---------------------------------------------------------------
 !Compute 1d spectra for various fields:
call spec1d(qqs,spec)
spec=log10(spmf*spec+1.d-32)
write(51,'(f12.5,1x,i5)') tgrid,kmaxred
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),spec(k)
enddo

call spec1d(zzs,spec)
spec=log10(spmf*spec+1.d-32)
write(52,'(f12.5,1x,i5)') tgrid,kmaxred
do k=1,kmaxred
  write(52,'(2(1x,f12.8))') alk(k),spec(k)
enddo

!---------------------------------------------------------------
 !Write scaled buoyancy (b_0/N) to bb.r4:
call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)
write(31,rec=igrids) real(tgrid),real(wka)

 !Write vertical vorticity to zz.r4:
call spctop(nx,ny,zzs,wka,xfactors,yfactors,xtrig,ytrig)
write(32,rec=igrids) real(tgrid),real(wka)

if (tracer) then
   !Write tracer (anomaly) spectrum cspec.asc:
  qqs=pt*cspre+ptc*cs
  call spec1d(qqs,spec)
  spec=log10(spmf*spec+1.d-32)
  write(53,'(f12.5,1x,i5)') tgrid,kmaxred
  do k=1,kmaxred
    write(53,'(2(1x,f12.8))') alk(k),spec(k)
  enddo

   !Write tracer (anomaly) field to cc.r4:
  call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)
  write(33,rec=igrids) real(tgrid),real(wka)
endif

write(*,'(a,f7.2,2(a,f13.6))') ' t = ',tgrid,' enstrophy = ',ens,' energy = ',ene

 !Unset flag for saving data:
gsave=.false.

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont

! Saves buoyancy contours for post-processing and imaging

implicit none

 !Local variables:
double precision:: ss(nx,ny)
double precision:: qa(ny,nx)
integer:: irec
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
irec=nint(t/tcsave)
write(pind(1:3),'(i3.3)') irec

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid buoyancy for plotting purposes:
ss=qq-qc
call spctop(nx,ny,ss,qa,xfactors,yfactors,xtrig,ytrig)
write(83,rec=irec) real(t),real(qa)

 !Save buoyancy contours if any exist:
if (contex) then
  open(81,file='cont/qqindex'//pind,form='unformatted',status='replace')
  write(81) npq(1:nq),i1q(1:nq),indq(1:nq)
  close(81)

  open(82,file='cont/qqnodes'//pind,form='unformatted',status='replace')
  write(82) xq(1:nptq),yq(1:nptq)
  close(82)
endif

 !Unset flag for saving data:
csave=.false.

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
