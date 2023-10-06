!#########################################################################
!  Finds balanced fields given only the PV anomaly in qq_init.r8 or
!  in qq.r8, using Optimal PV balance with fixed PV contours
!  (see Viudez & Dritschel, J. Fluid Mech. 521, 343-352, 2004).

!  The shape of the PV distribution is held fixed but its amplitude
!  is ramped up according to q(x,y,t)=qb(x,y)*R(t) where qb(x,y) is
!  the "base" configuration (read in from either qq_init.r8 or qq.r8),
!  while R(t) = 0.5*(1 - cos(pi*t/T)) is a ramp function going from
!  R = 0 at t = 0 to R = 1 at t = T, the specified ramping period.

!  This can be used both for initialisation, and for post-processing.

!   Adapted from VA code on 6 Feb 2020 by D G Dritschel @ St Andrews
!#########################################################################

program fopvbal

 ! Import constants and parameters:
use constants
 ! Import spectral module:
use spectral

implicit none
double precision:: hh(ng,ng),dd(ng,ng),gg(ng,ng),uu(ng,ng),vv(ng,ng)
double precision:: qb(ng,ng),qq(ng,ng),zz(ng,ng)
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)
double precision:: ekin,ediv,epot,etot,t,tramp
integer:: iopt,loop,iread,nstep

!---------------------------------------------------------
 ! Initialise inversion constants and arrays:
call init_spectral

write(*,*) ' Choose one of the following options:'
write(*,*)
write(*,*) ' (0) initialise from data in qq_init.r8'
write(*,*) ' (1) post-process from data in qq.r8'
write(*,*)
write(*,*) ' Option?'
read(*,*) iopt
write(*,*)
write(*,*) ' Ramping period, T?'
read(*,*) tramp
 ! Ensure tramp is a multiple of dt:
nstep=nint(tramp/dt)
tramp=dt*dble(nstep)

if (iopt .eq. 0) then
   ! Only find balance at t = 0 to initialise:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='old',recl=nbytes)
  read(11,rec=1) t,qb
  close(11)

   ! Find balanced fields:
  call balance

   ! Save data:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  write(11,rec=1) zero,qq
  close(11)

  open(11,file='dd_init.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  write(11,rec=1) zero,dd
  close(11)

  open(11,file='hh_init.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  write(11,rec=1) zero,hh
  close(11)

  open(11,file='gg_init.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  write(11,rec=1) zero,gg
  close(11)

  open(11,file='zz_init.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  write(11,rec=1) zero,zz
  close(11)

  write(*,*)
  write(*,*) ' The balanced fields are available in xx_init.r8 where'
  write(*,*) ' xx = dd, hh, gg, zz & qq (adjusted only by a constant).'
  write(*,*)
  
else
   ! Find balance at all times in qq.r8:
  open(30,file='evolution/qq.r8',form='unformatted', &
        access='direct',status='old',recl=nbytes)

  open(31,file='evolution/bqq.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(32,file='evolution/bdd.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(33,file='evolution/bhh.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(34,file='evolution/bgg.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(35,file='evolution/bzz.r8',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(51,file='evolution/becomp.asc',status='replace')

   ! Read data at all times and balance:
  loop=0
  do
    loop=loop+1
    iread=0
    read(30,rec=loop,iostat=iread) t,qb
    if (iread .ne. 0) exit

     ! Balance using only the PV anomaly at time t = t:
    write(*,'(a,f9.2)') ' *** Processing t = ',t
    call balance

     ! Write data:
    write(31,rec=loop) t,qq
    write(32,rec=loop) t,dd
    write(33,rec=loop) t,hh
    write(34,rec=loop) t,gg
    write(35,rec=loop) t,zz

    write(51,'(f9.2,5(1x,f16.9))') t,ediv,ekin,ekin+ediv,epot,etot
  enddo

  close(30)
  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(51)

  write(*,*)
  write(*,*) ' The balanced fields are available in bxx.r8 where xx = qq,'
  write(*,*) ' dd, hh, gg & zz.  Also, becomp.asc contains the energies.'
  write(*,*)

endif


 ! Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine balance

implicit none

 ! Local variables:
double precision:: htot(ng,ng)
integer:: istep

!-----------------------------------------------------------------
! Initialise:
ds=zero
gs=zero
qs=zero

dd=zero

!-----------------------------------------------------------------
! Perform the entire time integration:
do istep=1,nstep
  t=dt*dble(istep-1)
  write(*,'(a,f8.6)') ' t/T = ',t/tramp
   ! Evolve the flow from time t to t + dt:
  call advance
enddo
write(*,'(a,f8.6)') ' t/T = ',one

!-----------------------------------------------------------------
! Invert PV at final time:
call ptospc(ng,ng,qb,qs,xfactors,yfactors,xtrig,ytrig)
call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
call spctop(ng,ng,ds,dd,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,gs,gg,xfactors,yfactors,xtrig,ytrig)

! Compute energy components, divided by mean depth H:
htot=one+hh
ekin=f12*garea*sum(htot*(uu**2+vv**2))
epot=f12*garea*csq*sum(hh**2)
ediv=zero
etot=ekin+epot+ediv

return
end subroutine balance

!=======================================================================

subroutine advance

! Evolves delta & gamma_l by an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

! Local variables:
integer,parameter:: niter=2

! Spectral fields needed in time stepping:
double precision:: dsi(ng,ng),sds(ng,ng),nds(ng,ng)
double precision:: gsi(ng,ng),sgs(ng,ng),ngs(ng,ng)
! Other variables:
double precision:: ramp
integer:: iter

!-----------------------------------------------------------------------
! Invert PV, delta and gamma to obtain the velocity field:
call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
! Note: qs, ds & gs are in spectral space while 
!       hh, uu, vv, qq and zz are in physical space.

!------------------------------------------------------------------
! Start with a guess for F^{n+1} for all fields:

! Calculate the source terms (sds,sgs) for delta and h (ds,gs):
if (t .gt. zero) then
  call source(sds,sgs)
else
  sds=zero
  sgs=zero
endif

 !Update divergence and acceleration divergence:
dsi=ds
gsi=gs
nds=sds+dt4i*dsi
ngs=sgs+dt4i*gsi
sds=nds+sds          !2*N_tilde_delta
sgs=ngs+sgs          !2*N_tilde_gamma
ds=simp*(rdis*sds+sgs)-dsi       !simp = 1/(R^2-G);  rdis = R
gs=simp*(rdis*sgs+opak*sds)-gsi  !opak = G

! Update current PV at t+dt:
t=t+dt
ramp=ramp_fun(t)
qq=ramp*qb
call ptospc(ng,ng,qq,qs,xfactors,yfactors,xtrig,ytrig)

!------------------------------------------------------------------
! Iterate to improve estimates of F^{n+1}:
do iter=1,niter
  ! Perform inversion at t^{n+1} from estimated quantities:
  call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)

  ! Calculate the source terms:
  call source(sds,sgs)

   !Update divergence and acceleration divergence:
  sds=nds+sds          !2*N_tilde_delta
  sgs=ngs+sgs          !2*N_tilde_gamma
  ds=simp*(rdis*sds+sgs)-dsi       !simp = 1/(R^2-G);  rdis = R
  gs=simp*(rdis*sgs+opak*sds)-gsi  !pgop = G
enddo

return
end subroutine advance

!=======================================================================

subroutine source(sds,sgs)

! Gets the source terms (sds,sgs) for delta and gamma (ds,gs) ---
! all in spectral space.

! Note that (sds,sgs) only include the nonlinear terms for a 
! semi-implicit treatment, closely analogous to that described in 
! the appendix of Dritschel & Jalali (JFM, 2020).
  
implicit none

 ! Passed variables:
double precision:: sds(ng,ng),sgs(ng,ng)

 ! Local variables (physical):
double precision:: wkp(ng,ng),wkq(ng,ng)

 ! Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng)

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
end subroutine source

!=======================================================================

double precision function ramp_fun(time)

implicit none

double precision:: time

ramp_fun=f12-f12*cos(pi*time/tramp)

return
end function ramp_fun

 ! End main program
end program fopvbal
!=======================================================================
