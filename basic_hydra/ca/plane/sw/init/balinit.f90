!#########################################################################
!  Re-initialises a flow with balanced fields obtained from the conditions 
!  delta_t=gamma_t=0 using data previously set up with a data generation
!  routine.  Assumes the previous data has delta = gamma = 0.

!           Written 6/4/2018 by D G Dritschel @ St Andrews
!#########################################################################

program dgbalini

 !Import contants, parameters and common arrays:
use constants
use spectral

implicit none

 !Physical fields:
double precision:: qq(ng,ng),hh(ng,ng),dd(ng,ng),gg(ng,ng)
double precision:: uu(ng,ng),vv(ng,ng),zz(ng,ng),aa(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng),htot(ng,ng)
double precision:: hhpre(ng,ng),ddpre(ng,ng),ggpre(ng,ng),zzpre(ng,ng)

 !Spectral fields:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng),helm(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng),wke(ng,ng)

 !Other constants:
double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      hh, uu & vv from qq, dd & gg.  The energy error is computed from 
 !      <(u-u0)^2+(v-v0)^2+c^2*(h-h0)^2>/<u^2+v^2+c^2*h^2>
 !      where <:> means a domain average and (u0,v0,h0) is the previous
 !      guess for (u,v,h).
double precision:: qadd,qbar,fqbar,t,uio,vio
double precision:: ddrmserr,ggrmserr,toterr,toterrpre

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral
helm=-filt/opak

!----------------------------------------------------------------------
 !Read in gridded PV anomaly and convert to spectral space as qs:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)
 !Note: zz typically has zero domain average, whereas the actual
 !      PV anomaly may not since this is determined by the 
 !      requirement that the mean relative vorticity is zero;
 !      qr is corrected upon calling main_invert in spectral.f90

 !Convert to spectral space (zz is overwritten; the PV is recovered below):
call ptospc(ng,ng,zz,qs,xfactors,yfactors,xtrig,ytrig)

 !Ensure domain average qs is zero (this does not matter):
qs(1,1)=zero
 !Spectrally-truncate for use in de-aliasing:
qs=filt*qs

!----------------------------------------------------------------------
 !Start with zero divergence and acceleration divergence:
dd=zero
ds=zero
gg=zero
gs=zero

 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv),
 !relative vorticity (zz) and adjusted PV anomaly consistent with 
 !zero domain averaged zz:
call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
 !Note: qs, ds & gs are in spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

!-----------------------------------------------------------------
!Iterate to find the balanced fields:

!Energy norm error (must be > tole to start):
toterrpre=one/small**2
toterr=f12
hhpre=hh
ddpre=dd
ggpre=gg
zzpre=zz
do while (toterr .gt. tole)
  !Obtain balanced estimate for gamma (gg):
  call jacob(uu,vv,wkp)
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  wkp=dd*uu
  wkq=dd*vv
  call divs(wkp,wkq,wka)
  gs=filt*(wka-two*wkb)
  wka=gs
  call spctop(ng,ng,wka,gg,xfactors,yfactors,xtrig,ytrig)
  ggrmserr=sum((gg-ggpre)**2)/(sum(ggpre**2)+small)

  !Obtain balanced estimate for delta (dd):
  wkp=hh*uu
  wkq=hh*vv
  call divs(wkp,wkq,wka)
  wkp=zz*uu
  wkq=zz*vv
  call divs(wkp,wkq,wkb)
  ds=helm*(c2g2*wka-cof*wkb)
  wka=ds
  call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
  ddrmserr=sum((dd-ddpre)**2)/(sum(ddpre**2)+small)

  !Find height anomaly field (hh):
  htot=one+hh
  qadd=-dsumi*sum(qq*htot)
  qq=qq+qadd
  qbar=dsumi*sum(qq)
  wkp=cof*(qq+hh*(qq-qbar))-gg
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  fqbar=cof*qbar
  wka=filt*wkb/(opak-fqbar)
  call spctop(ng,ng,wka,hh,xfactors,yfactors,xtrig,ytrig)
   !wkp: corrected de-aliased height field (to be hh below)
  htot=one+hh

  !Obtain relative vorticity field (zz):
  wkp=qq*htot
  qadd=-dsumi*sum(wkp)
  qq=qq+qadd
  zz=htot*(cof+qq)-cof

  !Obtain velocity field (uu,vv):
  call ptospc(ng,ng,zz,wkb,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkb
  wkb=filt*wkb
  call spctop(ng,ng,wkb,zz,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)
  wke=rlap*ds
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkc)
  wkb=wka-wkb
  wkd=wkc+wkd
  call spctop(ng,ng,wkb,uu,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,vv,xfactors,yfactors,xtrig,ytrig)

  !Add mean flow (uio,vio):
  uio=-sum(hh*uu)*dsumi
  vio=-sum(hh*vv)*dsumi
  uu=uu+uio
  vv=vv+vio

  !Compute overall error:
  toterr=f12*(ddrmserr+ggrmserr)

  write(*,*) ' relative delta error = ',ddrmserr
  write(*,*) ' relative gamma error = ',ggrmserr

  !If error is going up again, stop and save fields:
  if (toterrpre .lt. toterr) exit

  !Otherwise continue with another iteration:
  hhpre=hh
  ddpre=dd
  ggpre=gg
  zzpre=zz
  toterrpre=toterr
enddo

write(*,*) ' Minimum error = ',toterrpre

!-----------------------------------------------------------------
!Write data:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,dd
close(11)

open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,gg
close(11)

write(*,*)
write(*,*) ' Initial fields balanced and re-written.'

 !End main program
end program
!=======================================================================
