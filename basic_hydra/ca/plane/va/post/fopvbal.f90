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

!   Adapted from 2-layer code 1 Feb 2020 by D G Dritschel @ St Andrews
!#########################################################################

program fopvbal

 ! Import constants and parameters:
use constants
 ! Import spectral module:
use spectral

implicit none
double precision:: hh(ng,ng),dd(ng,ng),gg(ng,ng),uu(ng,ng),vv(ng,ng)
double precision:: qb(ng,ng),qq(ng,ng),zz(ng,ng),ppn(ng,ng)
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

  open(36,file='evolution/bpn.r8',form='unformatted', &
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
    write(36,rec=loop) t,ppn

    write(51,'(f9.2,5(1x,f16.9))') t,ediv,ekin,ekin+ediv,epot,etot
  enddo

  close(30)
  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(36)
  close(51)

  write(*,*)
  write(*,*) ' The balanced fields are available in bxx.r8 where xx = qq,'
  write(*,*) ' dd, hh, gg, zz & ppn.  Also, becomp.asc contains the energies.'
  write(*,*)

endif


 ! Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine balance

implicit none

 ! Local variables:
double precision:: htot(ng,ng)
double precision:: wka(ng,ng)
integer:: istep

!-----------------------------------------------------------------
! Initialise:
ds=zero
gs=zero
qs=zero

ppn=zero
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
wka=htot*dd
ediv=f12*garea*hbsq3*sum(htot*wka**2)
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

! Update divergence and acceleration divergence:
dsi=ds
gsi=gs
nds=sds+dt4i*dsi
ngs=sgs+dt4i*gsi
sds=nds+sds                      !2*N_tilde_delta
sgs=ngs+sgs                      !2*N_tilde_gamma
ds=simp*(pdis*sds+sgs)-dsi       !simp = 1/(P*R^2-G);  pdis = P*R
gs=simp*(pdis*sgs+pgop*sds)-gsi  !pgop = P*G

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

  ! Update divergence and acceleration divergence:
  sds=nds+sds                      !2*N_tilde_delta
  sgs=ngs+sgs                      !2*N_tilde_gamma
  ds=simp*(pdis*sds+sgs)-dsi       !simp = 1/(P*R^2-G);  pdis = P*R
  gs=simp*(pdis*sgs+pgop*sds)-gsi  !pgop = P*G
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

double precision,parameter:: toler=1.d-11
 ! toler: maximum error in iteration below to find the vertically-
 !        integrated non-hydrostatic pressure \bar{P}_n (ppn below).
double precision,parameter:: w=0.75d0, wc=one-w
 ! w: weight used in the pressure iteration to accelerate convergence

 ! Passed variables:
double precision:: sds(ng,ng),sgs(ng,ng)

 ! Local variables (physical):
double precision:: htot(ng,ng),hinv(ng,ng),hinvsq(ng,ng)
double precision:: rx(ng,ng),ry(ng,ng),pnx(ng,ng),pny(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng)
double precision:: errpn

 ! Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng),wkf(ng,ng),wkg(ng,ng)

!---------------------------------------------------------------
 ! Calculate the vertically-integrated non-hydrostatic pressure.

 ! First form gamma_tilde = gamma + 2*(J(u,v) - delta^2):
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
 ! dd contains the divergence (delta) in physical space

call jacob(uu,vv,wkp)
 ! wkp contains J(u,v) in physical space

wkp=wkp-dd**2
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
wkg=gs+two*filt*wkg
call spctop(ng,ng,wkg,wkp,xfactors,yfactors,xtrig,ytrig)
 ! wkp now contains gamma_tilde in physical space (de-aliased)

 ! Multiply next by (1+h) and re-store in wkg (spectral) as the 
 ! fixed rhs in the pressure iteration below:
htot=one+hh
wkp=htot*wkp
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
 ! wkg is not de-aliased by the prop operator below takes care of this.

 ! Next calculate wkq = 3/H^2*(1/(1+h)^2 - 1) needed for the pressure
 ! iteration below:
hinv=one/htot
call dealias(hinv)
hinvsq=hinv**2
call dealias(hinvsq)
wkq=hbsq3i*(hinvsq-one)

 ! Calculate also (1+h)^{-1}*grad(h) and store in rx & ry:
wkp=hh
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,rx,ry)
rx=hinv*rx
call dealias(rx)
ry=hinv*ry
call dealias(ry)

 ! Now iterate to find \bar{P}_n (in ppn) starting from the guess
 ! ppn = (grad^2 - 3/H^2)^{-1}((1+h)*gamma_tilde):
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

 ! Store rhs of pressure iteration, minus gs, in wkb for use in divergence
 ! tendency below:
wkb=wkg-gs+wkb

 ! Compute J(ppn,hinv) and store spectral version in wkf (not de-aliased):
call jacob(ppn,hinv,wkp)
call ptospc(ng,ng,wkp,wkf,xfactors,yfactors,xtrig,ytrig)
 ! wkf is needed below for the gamma source.

!---------------------------------------------------------------
 ! Nonlinear part of ds source

 ! Compute div(delta*u,delta*v) and store in wka (spectral):
wkp=dd*uu
wkq=dd*vv
call divs(wkp,wkq,wka)

 ! Compute delta^2 and store in spectral space as wkg:
wkp=dd**2
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)

 ! Store (spectral) 2*delta^2 - div(delta*u,delta*v) in wka to free up wkg:
wka=two*wkg-wka

 ! Compute 1/(1+h)^3 -> wkp and de-alias:
wkp=hinv*hinvsq
call dealias(wkp)
 ! Form \bar{P}_n*(1/(1+h)^3 - 1):
wkp=ppn*(wkp-one)
 ! Transform to spectral space as wkg:
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)

 ! Add everything up to define delta_t - pope*gamma (spectral, filtered)
sds=filt*(wka-hbsq3i*wkg+pope*wkb)
 ! Here pope = (1 - (H^2/3)*grad^2)^{-1} (spectral)

!---------------------------------------------------------------
 ! Nonlinear part of gs source --- f*N_zeta - c^2*Lap(N_h) where
 ! N_zeta and N_h are the nonlinear parts of the vorticity and
 ! depth sources:

 ! Obtain N_h:
wkp=hh*uu
wkq=hh*vv
 ! Compute div(h*u,h*v) = -N_h spectrally and filter for use below:
call divs(wkp,wkq,wka)
 ! For use below, calculate -c^2*Lap(N_h):
wka=c2g2*wka

 ! Obtain N_zeta:
wkp=zz*uu
wkq=zz*vv
 ! Compute div(zeta*u,zeta*v) spectrally:
call divs(wkp,wkq,wkb)

 ! Add everything up to define S_gamma:
sgs=wka+cof*filt*(wkf-wkb)
 ! wkf contains J(h,2*P_0/3)/(1+h) in spectral space, part of N_zeta.

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
