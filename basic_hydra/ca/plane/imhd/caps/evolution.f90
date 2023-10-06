module evolution

! Module contains subroutines to evolve PV contours and all fields 
! according to the algorithm detailed in caps.f90.

use common

implicit none

 !Velocity field:
double precision:: uu(ny,nx),vv(ny,nx)

 !Spectral fields (note array order):
double precision:: pp(nx,ny),qq(nx,ny),qc(nx,ny),qd(nx,ny)
double precision:: qspre(nx,ny),aapre(nx,ny)
double precision:: emq(nx,ny),epq(nx,ny),ema(nx,ny),epa(nx,ny)

 !Energies:
double precision:: eneupre,enebpre

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
double precision,parameter:: qratmax=0.2d0
!      qratmax:  the maximum ratio r of the mean-square residual PV qd
!                to the mean-square PV qc in the contours
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, or when r > qratmax, 
!      the code rebuild the PV contours in a separate memory space.
double precision:: qa(ny,nx),qrat
double precision:: wka(nx,ny)
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
    call spctop(nx,ny,wka,qr,xfactors,yfactors,xtrig,ytrig)
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
    write(14,'(1x,f12.5,1x,i8,1x,i9)') t,na,npta

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!Save final data if not done already:
if (gsave) call savegrid
if (csave) call savecont

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises quantities needed for normal time integration following 
! contour regeneration

implicit none

 !Local variables:
double precision:: qa(ny,nx)
double precision:: wka(nx,ny),qmin,ql1,ql2
integer:: nc,nptc,i,j,ix,iy

!------------------------------------------------------------------
 !Define number of "active" contours and nodes:
na=nq
npta=nptq

if (tracer) then
   !There are tracer contours; read them from tracer.bin and append
   !them on to the list of contours:
  open(11,file='contours/tracer.bin',status='old',form='unformatted')
  read(11) nc,nptc

  nq=na+nc
  do j=na+1,nq
    read(11) npq(j)
    i1q(j)=nptq+1
    i2q(j)=nptq+npq(j)
    nptq=i2q(j)
    indq(j)=9999
  enddo
   !Tracer contours are assigned an index of 9999 to distinguish them
   !from active contours.

  do i=npta+1,nptq
    read(11) xq(i),yq(i)
    nextq(i)=i+1
  enddo

  close(11)

   !Correct nextq array, giving the node following a node:
  do j=na+1,nq
    nextq(i2q(j))=i1q(j)
  enddo
endif

 !Logical to indicate presence of contours:
contex=(nptq .gt. 0)

 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,na,npta

!-------------------------------------------------------------
 !Logicals used for saving gridded fields and contours:
gsave=.false.
csave=.false.

 !Convert PV contours to gridded values (qa):
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define residual PV qd = (1-F)[qs-qc], and copy current gridded 
 !fields (qs,aa) for use in time interpolation:
qd=fhi*(qs-qc)
qs=filt*qs 
qspre=qs 
aa=filt*aa 
aapre=aa
 !Here fhi = 1-F is a high-pass spectral filter and filt is a 
 !de-aliasing filter (defined in spectral.f90)

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current PV anomaly field is stored in qs, the residual PV needed 
! in congen.f90 is stored in qr, and (if present) tracer contours are
! stored in tracer.bin

implicit none

 !Local variables:
double precision:: qa(ny,nx)
double precision:: wka(nx,ny),qmin,ql1,ql2
integer:: i,j,ix,iy

!-----------------------------------------------------------------
 !Convert PV contours to gridded values (qa):
call con2grid(qa)

 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Put current PV anomaly (spectral in qs) and define residual (qd):
qs=flo*qs+fhi*qc+qd
qd=qs-qc

 !Convert qd to physical space as qr (used in recontouring):
call spctop(nx,ny,qd,qr,xfactors,yfactors,xtrig,ytrig)

!-------------------------------------------------------------------
 !Update PV contour interval for re-contouring:
if (beffect) then
   !Choose an integral number of PV jumps (required in recontour):
  qjump=beta*hly/dble(ncontq)
else
   !Choose contour interval based on <q^2>/<|q|> for |q| > q_rms:
  wka=qs
  call spctop(nx,ny,wka,qa,xfactors,yfactors,xtrig,ytrig)
  qmin=sqrt(dsumi*sum(qa**2))
  ql1=zero
  ql2=zero
  do ix=1,nx
    do iy=1,ny
      if (abs(qa(iy,ix)) .gt. qmin) then
        ql1=ql1+abs(qa(iy,ix))
        ql2=ql2+qa(iy,ix)**2
      endif
    enddo
  enddo
  qjump=ql2/(ql1*dble(ncontq))
endif

!-------------------------------------------------------------------
if (tracer) then
   !There are tracer contours; write them to tracer.bin to be used
   !again when this module is re-entered after contour regeneration:
  open(11,file='contours/tracer.bin',status='replace',form='unformatted')
  write(11) nq-na,nptq-npta
  do j=na+1,nq
    write(11) npq(j)
  enddo
  do i=npta+1,nptq
    write(11) xq(i),yq(i)
  enddo
  close(11)
   !Only build new contours using active contours:
  nq=na
  nptq=npta
endif

return
end subroutine prepare

!=======================================================================

subroutine advance

! Advances PV and A from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method (for all
! gridded fields, i.e. qs, qd & aa).

! *** Uses a 4th-order Runge-Kutta method ***

implicit none

 !Local variables:

 !Spectral fields needed in Runge-Kutta time stepping (note array order):
double precision:: qsi(nx,ny),qsf(nx,ny),sqs(nx,ny)
double precision:: qdi(nx,ny),qdf(nx,ny),sqd(nx,ny)
double precision:: aai(nx,ny),aaf(nx,ny),saa(nx,ny)
 !Contour positions needed in Runge-Kutta time stepping:
double precision:: xqi(nptq),yqi(nptq),xqf(nptq),yqf(nptq)
 !Contour velocities:
double precision:: uq(nptq),vq(nptq)
 !Other local quantities:
double precision:: xx,yy
integer:: i

!-------------------------------------------------------------------
 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*qs + (1-F)*qc + qd
 !            and qd = (1-F)*(qs-qc)
 !after qs is reset; here F is a low pass filter (see spectral.f90)
qs=flo*qs+fhi*qc+qd
qd=fhi*(qs-qc)

!------------------------------------------------------------------
 !RK4 predictor step to time t0 + dt/2:

 !Invert PV and compute velocity:
call inversion

 !Possibly save data (gsave & csave set by adapt in the previous time step):
if (gsave) call savegrid
if (csave) call savecont

 !Adjust timestep (dt) on maximum vorticity magnitude or CFL:
call adapt

 !Calculate the source terms (sqs,sqd,saa) for PV (qs,qd) and A (aa):
call source(sqs,sqd,saa,0)

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
aai=aa
aa=ema*(aai+dt2*saa)
aaf=aai+dt6*saa

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt/2:
t=t+dt2

 !Invert PV and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd,saa) for PV (qs,qd) and A (aa):
call source(sqs,sqd,saa,1)

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
aa=ema*(aai+dt2*saa)
aaf=aaf+dt3*saa

!------------------------------------------------------------------
 !RK4 predictor step at time t0 + dt:

 !Invert PV and compute velocity:
call inversion

 !Update any vorticity and magnetic potential stochastic forcing:
if (zforcing) call zforce
if (aforcing) call aforce
 !Note, the 1st two substeps use dzdt & dadt at time t0,
 !while the 2nd two substeps use dzdt & dadt at time t0 + dt.

 !Calculate the source terms (sqs,sqd,saa) for PV (qs,qd) and A (aa):
call source(sqs,sqd,saa,1)

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
ema=ema**2
aa=ema*(aai+dt*saa)
aaf=aaf+dt3*saa

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt:
t=t+dt2

 !Invert PV and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd,saa) for PV (qs,qd) and A (aa):
call source(sqs,sqd,saa,2)

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
aa=ema*(aaf+dt6*saa)

return
end subroutine advance

!=======================================================================

subroutine source(sqs,sqd,saa,lev)

! Gets the source terms (sqs,sqd,saa) for the PV (qs,qd) & A (aa) 
! evolution equations (all in spectral space):

implicit none

 !Passed variables:
double precision:: sqs(nx,ny),sqd(nx,ny),saa(nx,ny)
integer:: lev

 !Local variables:
double precision:: aax(ny,nx),aay(ny,nx),qqx(ny,nx),qqy(ny,nx)
double precision:: wkp(ny,nx)
integer:: ix,iy

!---------------------------------------------------------------
 !qd source - first compute curl of Lorentz force (wkp):
call lorentz(aa,aax,aay,wkp)

 !Add on NL term for qd:
call gradient(qd,qqx,qqy)
wkp=wkp-uu*qqx-vv*qqy

 !Add stochastic vorticity forcing (if present):
if (zforcing) wkp=wkp+dzdt

 !Convert wkp to spectral space as sqd:
call ptospc(nx,ny,wkp,sqd,xfactors,yfactors,xtrig,ytrig)

 !Implement Ekman and/or thermal damping (if present):
if (heating) then
  if (friction) then 
   !Use thermal and Ekman damping:
    sqd=sqd+therm*pp-rekman*(qq+kdsq*pp)
  else
    sqd=sqd+therm*pp
  endif
else
   !Only use Ekman damping
  sqd=sqd-rekman*(qq+kdsq*pp)
endif

!---------------------------------------------------------------
 !qs source - only NL term is needed:
call gradient(qs,qqx,qqy)
wkp=-uu*qqx-vv*qqy

 !Convert to spectral space:
call ptospc(nx,ny,wkp,sqs,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !aa source - only NL term is needed:
wkp=-uu*aax-vv*(b0+aay)

 !Add stochastic magnetic potential forcing (if present):
if (aforcing) wkp=wkp+dadt

 !Convert wkp to spectral space as saa:
call ptospc(nx,ny,wkp,saa,xfactors,yfactors,xtrig,ytrig)

!----------------------------------------------------------------
if (lev .eq. 0) then
 !Spectrally truncate sources:
  sqs=filt*sqs
  sqd=filt*sqd
  saa=filt*saa
 !Apply exponential integrating factors (and spectrally truncate sqs & sqd):
else if (lev .eq. 1) then
  sqs=filt*sqs
  sqd= epq*sqd
  saa= epa*saa
else
  sqs=  filt*sqs
  sqd=epq**2*sqd
  saa=epa**2*saa
endif

return
end subroutine source

!=======================================================================

subroutine zforce

! Updates stochastic vorticity forcing dzdt (time correlated)

implicit none

 !Local variables:
double precision:: wkp(ny,nx)
double precision:: fac

!------------------------------------------------------------------
 !Generate a random field with a particular spectrum and rms value:
call ranspec(wkp,esrz,powz,k0z)

 !Blend with current forcing to get updated forcing at time t + dt:
fac=dt/tcz
dzdt=dzdt+fac*(wkp-dzdt)

 !Restore rms value of blended field:
fac=esrz/sqrt(dsumi*sum(dzdt**2))     !dsumi = 1/dble(nx*ny)
dzdt=fac*dzdt

return
end subroutine zforce

!=======================================================================

subroutine aforce

! Updates stochastic magnetic potential forcing dadt (time correlated)

implicit none

 !Local variables:
double precision:: wkp(ny,nx)
double precision:: fac

!------------------------------------------------------------------
 !Generate a random field with a particular spectrum and rms value:
call ranspec(wkp,esra,powa,k0a)

 !Blend with current forcing to get updated forcing at time t + dt:
fac=dt/tca
dadt=dadt+fac*(wkp-dadt)

 !Restore rms value of blended field:
fac=esra/sqrt(dsumi*sum(dadt**2))     !dsumi = 1/dble(nx*ny)
dadt=fac*dadt

return
end subroutine aforce

!=======================================================================

subroutine inversion

! Inverts Laplace's operator on PV anomaly (PV - beta*y) to obtain 
! the streamfunction (pp) ***in spectral space*** and the velocity 
! (uu,vv) = (-dpp/dy,dpp/dx) ***in physical space***

implicit none

 !Local variables:
double precision:: qa(ny,nx)

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Combine fields to update qq with full field,
 !qq = F*qs + (1-F)*qc + qd, where F is a low pass filter:
qq=flo*qs+fhi*qc+qd

 !Invert PV to obtain velocity field: 
call main_invert(qq,uu,vv,pp)

return
end subroutine inversion

!=======================================================================

subroutine adapt

! Adapts the time step dt to ensure that it is less than or equal to
! the minimum of dtfac/max(|zeta|_max) and C*dx/max(|u|_max,|b|_max)
! where dx is the grid spacing, u is the vector velocity field, and
! b is the vector magnetic field.  C = cfl_max is specified below.

implicit none

 !Local variables:
double precision,parameter:: dtfac=pi/10.d0, cflmax=0.7d0
double precision:: aax(ny,nx),aay(ny,nx),zz(ny,nx),ff(ny,nx) !Physical
double precision:: ss(nx,ny)                                 !Spectral
double precision:: umax,bmax,zzrms,zzmax,dtacc,tcont,epot,dfac
double precision:: circh,circm,dchdt
integer:: itime

!----------------------------------------------------------------------
 !Compute the gridded relative vorticity (zz):
if (barotropic) then
   !Here kd = 0:
  ss=qq
else
   !Here kd > 0:
  ss=qq+kdsq*pp
endif
 !Above, qq = PV anomaly (q - beta*y) in spectral space.
call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)

 !Compute the gradients of A to get B = (B_0 + A_y, -A_x)
call gradient(aa,aax,aay)

 !Compute accurate advection time step:
umax=sqrt(maxval(uu**2+vv**2))+small
bmax=sqrt(maxval(aax**2+(b0+aay)**2))
ff=zz**2
zzmax=sqrt(maxval(ff))+small
zzrms=sqrt(dsumi*sum(ff))
dtacc=min(glx*cflmax/max(umax,bmax),dtfac/max(zzmax,srwfm))
 !The restriction on the maximum Rossby wave frequency (srwfm)
 !ensures that the fastest Rossby wave frequency is resolved.

 !Write rms vorticity, max vorticity, max|u| & |B| to monitor.asc:
write(17,'(1x,f12.5,1x,f12.7,1x,f14.7,2(1x,f12.7))') t,zzrms,zzmax,umax,bmax

!---------------------------------------------------------------------
 !Choose a new time step, dt:
dt=min(dtacc,dtmax)
if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc

 !Fractional time steps used in 4th-order Runge-Kutta time stepping:
dt2=dt*f12
dt3=dt*f13
dt6=dt*f16

 !Increment the integral of max|zz|:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
if (itime+1 .gt. igrids) then
   !The save time is between t & t+dt; set flag to save data:
  gsave=.true.

   !The next gridded data save time: 
  tgrid=tgsave*dble(itime)

   !Copy current gridded fields for use in time interpolation:
  qspre=qs 
  aapre=aa

   !Compute fluid and magnetic energy:
  call energy(aax,aay,eneupre,enebpre)
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
emq=exp(-dfac*qdiss)
epq=filt/emq
ema=exp(-dt2*adiss)
epa=filt/ema

!----------------------------------------------------------------------
 !Compute circulation statistics if tracer contours are present:
if (tracer) then
   !Define j in spectral space:
  ss=rksq*aa
  call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)
   !Use the array zz to hold the gridded current density
  call circulation(uu,vv,aax,aay,zz,circh,circm,dchdt)
   !Record results in circul.dat:
  write(24,'(f12.5,1x,f14.9,2(1x,f17.12))') t,circh,circm,dchdt
endif

return
end subroutine adapt

!=======================================================================

subroutine energy(aax,aay,eu,eb)

! Given (A_x,A_y) in (aax,aay), this routine computes the total
! fluid and magnetic energies from uu, vv & pp (if kd > 0).

implicit none

 !Passed variables:
double precision:: aax(ny,nx),aay(ny,nx) !Physical
double precision:: eu,eb

 !Local variables:
double precision:: ss(nx,ny) !Spectral
double precision:: zz(nx,ny) !Physical
double precision:: epot

!----------------------------------------------------------------------
if (barotropic) then
   !kd = 0 and hence there is no potential energy:
  epot=zero
else
   !Get streamfunction pp in physical space as zz:
  ss=pp
  call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)
   !Compute potential energy per unit area:
  epot=f12*dsumi*kdsq*sum(zz**2)
endif

 !Compute kinetic and magnetic energies per unit area:
zz=uu**2+vv**2
eu=f12*dsumi*sum(zz)+epot
zz=aax**2+aay**2
eb=f12*dsumi*sum(zz)

return
end subroutine energy

!=======================================================================

subroutine savegrid

! Saves PV, j, A, energy and various spectra at the desired save time

implicit none

 !Local variables:
double precision:: wka(ny,nx),wkb(ny,nx) !Physical
double precision:: aas(nx,ny),qqs(nx,ny) !Spectral
double precision:: qspec(0:max(nx,ny)),jspec(0:max(nx,ny))
double precision:: qql2,jjl2,aal2,pt,ptc,eneu,eneb
double precision:: eneupost,enebpost
real:: tr4
integer:: ix,k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Interpolate spectral magnetic potential at save time:
aas=pt*aapre+ptc*aa

 !Write gridded A to a file:
qqs=aas   !use qqs temporarily to avoid overwriting aas
call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)
tr4=real(tgrid)
write(33,rec=igrids) tr4,real(wka)

 !Compute mean-square A:
aal2=dsumi*sum(wka**2)

 !Interpolate spectral PV anomaly at save time:
qqs=pt*qspre+ptc*qq

 !Compute the gradient of A and put it in (wka,wkb):
call gradient(aa,wka,wkb)

 !Compute fluid and magnetic energy:
call energy(wka,wkb,eneupost,enebpost)

 !Compute time interpolated energies:
eneu=pt*eneupre+ptc*eneupost
eneb=pt*enebpre+ptc*enebpost

 !Write energy to ene.asc:
write(15,'(f12.5,3(1x,f16.9))') tgrid,eneu,eneb,eneu+eneb

!---------------------------------------------------------------
 !Compute 1d PV and current density spectra:
call spec1d(qqs,qspec,0)
call spec1d(aas,jspec,1)
 !*** Note: aas is returned as j in spectral space

 !Write PV spectrum:
write(51,'(f12.5,1x,i5)') tgrid,kmaxred
qspec=log10(spmf*qspec+1.d-32)
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),qspec(k)
enddo

 !Write current density spectrum:
write(52,'(f12.5,1x,i5)') tgrid,kmaxred
jspec=log10(spmf*jspec+1.d-32)
do k=1,kmaxred
  write(52,'(2(1x,f12.8))') alk(k),jspec(k)
enddo

!----------------------------------------------------------------------
 !Compute domain average of (q-beta*y)^2 & j^2 and write to norms.asc:
call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)
qql2=dsumi*sum(wka**2)

call spctop(nx,ny,aas,wkb,xfactors,yfactors,xtrig,ytrig)
jjl2=dsumi*sum(wkb**2)
 !Note: dsumi = 1/(nx*ny)

write(16,'(f12.5,3(1x,f16.9))') tgrid,qql2,jjl2,aal2

write(*,'(a,f7.2,2(a,f9.4),2(a,f8.6))') ' t = ',tgrid, &
   & ' <q^2> = ',qql2,' <j^2> = ',jjl2,' E_u = ',eneu,' E_b = ',eneb

!----------------------------------------------------------------------
 !Write full PV to qq.r4 and current density to jj.r4:
if (beffect) then
   !Add beta*y to define full PV:
  do ix=1,nx
    wka(:,ix)=wka(:,ix)+bety
  enddo
endif
write(31,rec=igrids) tr4,real(wka)
write(32,rec=igrids) tr4,real(wkb)

 !Unset flag for saving data:
gsave=.false.

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont

! Saves PV contours for post-processing and imaging

implicit none

 !Local variables:
double precision:: ss(nx,ny)
double precision:: qa(ny,nx)
integer:: irec
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f7.2)') ' Saving contours at t = ',t
irec=nint(t/tcsave)
write(pind(1:3),'(i3.3)') irec

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
ss=qq-qc
call spctop(nx,ny,ss,qa,xfactors,yfactors,xtrig,ytrig)
write(83,rec=irec) real(t),real(qa)

 !Save PV contours if any exist:
if (contex) then
  open(81,file='contours/qqindex'//pind,form='unformatted', &
        access='direct',status='replace',recl=12*nq)
  write(81,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq)
  close(81)

  open(82,file='contours/qqnodes'//pind,form='unformatted', &
        access='direct',status='replace',recl=16*nptq)
  write(82,rec=1) xq(1:nptq),yq(1:nptq)
  close(82)
endif

 !Unset flag for saving data:
csave=.false.

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
