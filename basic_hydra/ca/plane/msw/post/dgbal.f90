!######################################################################
!  Finds balanced fields given only the PV anomaly in qq_init.r8 or
!  in qq.r8, using first-order delta-gamma balance in which delta_t
!  and gamma_t are set to zero to diagnose all fields.  Here gamma
!  is the linearised acceleration divergence.

!  This can be used both for initialisation and for post-processing.

!  Adapted from VA code on 6 Feb 2020 by D G Dritschel @ St Andrews
!######################################################################

program dgbal

 ! Import constants and parameters:
use constants
 ! Import spectral module:
use spectral

implicit none
double precision:: hh(ng,ng),dd(ng,ng),gg(ng,ng),uu(ng,ng),vv(ng,ng)
double precision:: qb(ng,ng),qq(ng,ng),zz(ng,ng)
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)
double precision:: ekin,ediv,epot,etot,t,berro
integer:: iopt,loop,iread,nstep

!---------------------------------------------------------
 ! Initialise inversion constants and arrays:
call init_spectral

 ! Redefine operators needed to solve for the balanced fields:
simp=-filt/opak

write(*,*) ' Choose one of the following options:'
write(*,*)
write(*,*) ' (0) initialise from data in qq_init.r8'
write(*,*) ' (1) post-process from data in qq.r8'
write(*,*)
write(*,*) ' Option?'
read(*,*) iopt
write(*,*)

if (iopt .eq. 0) then
   ! Only find balance at t = 0 to initialise:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='old',recl=nbytes)
  read(11,rec=1) t,qb
  close(11)

   ! Find balanced fields:
  call balance(berro)

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

  open(61,file='evolution/benorm.asc',status='replace')

   ! Read data at all times and balance:
  loop=0
  do
    loop=loop+1
    iread=0
    read(30,rec=loop,iostat=iread) t,qb
    if (iread .ne. 0) exit

     ! Balance using only the PV anomaly at time t = t:
    write(*,'(a,f9.2)') ' *** Processing t = ',t
    call balance(berro)

     ! Write data:
    write(31,rec=loop) t,qq
    write(32,rec=loop) t,dd
    write(33,rec=loop) t,hh
    write(34,rec=loop) t,gg
    write(35,rec=loop) t,zz

    write(51,'(f9.2,5(1x,f16.9))') t,ediv,ekin,ekin+ediv,epot,etot

    write(61,'(f9.2,1x,f11.7)') t,log10(berro)
  enddo

  close(30)
  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(51)
  close(61)

  write(*,*)
  write(*,*) ' The balanced fields are available in bxx.r8 where xx = qq,'
  write(*,*) ' dd, hh, gg & zz.'
  write(*,*)
  write(*,*) ' Also, becomp.asc contains the energy components, while'
  write(*,*) ' benorm.asc contains log10(relative energy norm error).'
  write(*,*)

endif


 ! Internal subroutine definitions (inherit global variables):

contains

!=======================================================================

subroutine balance(berro)

implicit none

! Passed variable (returned):
double precision:: berro

! Local variable:
double precision:: htot(ng,ng)

!-----------------------------------------------------------------
! Iterate to find balanced fields:
call iterate(berro)

!-----------------------------------------------------------------
! Invert PV again after convergence:
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

subroutine iterate(berro)

! Solves the equations resulting from setting delta_t = gamma_t = 0.
!
! Exploits the semi-implicit formulation, here modified for an
! infinite time step and no dissipation.

implicit none

! Error in energy norm below which solution is said to be converged:
double precision,parameter:: tole=1.d-14

! Passed variable (returned):
double precision:: berro

! Local variables:
double precision:: hho(ng,ng),uuo(ng,ng),vvo(ng,ng)
double precision:: dso(ng,ng),gso(ng,ng)
double precision:: sds(ng,ng),sgs(ng,ng)
double precision:: berr

!-----------------------------------------------------------------------
! Initialise:
ds=zero
gs=zero

dso=zero
gso=zero

hho=zero
uuo=zero
vvo=zero
qq=qb
call ptospc(ng,ng,qq,qs,xfactors,yfactors,xtrig,ytrig)

! Iterate to find balance (minimise error between two iterates):
berro=one
do
   ! Invert PV, delta and gamma to obtain the velocity field:
  call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
   ! Note: qs, ds & gs are in spectral space while 
   !       hh, uu, vv, qq and zz are in physical space.

   ! Compute energy norm error:
  berr=sum((uu-uuo)**2+(vv-vvo)**2+csq*(hh-hho)**2)/sum(uu**2+vv**2+csq*hh**2)

  write(*,*) ' Energy error = ',berr

  if (berr .lt. tole) then
     ! Error below minimum specified:
    berro=berr
    write(*,*) ' Converged!  Relative energy error = ',berro
     ! Exit the loop and return from this subroutine
    exit
  endif

  if (berr .gt. berro) then
     ! Error increasing again; use previous iterate as the balanced state:
    ds=dso
    gs=gso
    write(*,*) ' Minimum relative energy error = ',berro
     ! Exit the loop and return from this subroutine
    exit
  endif

  ! Copy new fields into old ones for next iteration:
  uuo=uu
  vvo=vv
  hho=hh
  dso=ds
  gso=gs
  berro=berr

  ! Compute NL source terms for delta and gamma:
  call source(sds,sgs)

   ! Update divergence and acceleration divergence:
  ds=simp*sgs
  gs=-sds
enddo

return
end subroutine iterate

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

 ! End main program
end program dgbal
!=======================================================================
