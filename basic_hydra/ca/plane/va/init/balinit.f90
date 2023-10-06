!#########################################################################
!  Re-initialises a flow with balanced fields obtained from the conditions 
!  delta_t=gamma_t=0 using data previously set up with a data generation
!  routine.  Assumes the previous data has delta = gamma = 0.

!           Written 6/4/2018 by D G Dritschel @ St Andrews
!#########################################################################

program dgbalini

 !Import contants, parameters and common arrays:
use spectral

implicit none

double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      hh, uu & vv from qq, dd & gg.  The energy error is computed from 
 !      <(u-u0)^2+(v-v0)^2+c^2*(h-h0)^2>/<u^2+v^2+c^2*h^2>
 !      where <:> means a domain average and (u0,v0,h0) is the previous
 !      guess for (u,v,h).

 !Physical fields:
double precision:: qq(ng,ng),hh(ng,ng),dd(ng,ng),gg(ng,ng)
double precision:: uu(ng,ng),vv(ng,ng),zz(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng)
double precision:: htot(ng,ng),hx(ng,ng),hy(ng,ng)
double precision:: hhpre(ng,ng),uupre(ng,ng),vvpre(ng,ng)

 !Spectral fields:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng),wke(ng,ng)

 !Other constants:
double precision:: uio,vio,t
double precision:: qadd,qbar,fqbar
double precision:: dhrms,durms,enorm

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!----------------------------------------------------------------------
 !Read in approximate PV anomaly, (zeta+f)/(1+h_tilde) and convert to 
 !spectral space as qs:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
read(11,rec=1) t,qq
close(11)

 !Convert to spectral space (zz is overwritten):
zz=qq
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

 !Initialise previous guess for hh, uu & vv:
hh=zero
uu=zero
vv=zero

 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv),
 !relative vorticity (zz) and adjusted PV anomaly consistent with 
 !zero domain averaged zz:

 !Define first guess for total dimensionless height:
htot=one
!-------------------------------------------------------
 !Iteratively solve for hh, uu & vv:

 !Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)
   !Get average PV from the requirement of zero average vorticity:
  qadd=-dsumi*sum(qq*htot)
   !qadd+qq is the corrected PV anomaly; dsumi = 1/ng^2 here.

  qq=qq+qadd
  qbar=dsumi*sum(qq)
   !qbar = <q> is the average PV anomaly; use to invert for h:

   !Invert [c^2*grad^2-f*(f+<q>)]h = f*(q+h*(q-<q>):
  wkp=cof*(qq+hh*(qq-qbar))
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)

  fqbar=cof*qbar
  wka=filt*wkb/(opak-fqbar)
   !opak: c^2*grad^2-f^2 (spectral)
  call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
   !wkp: corrected de-aliased height field (to be hh below)

   !Compute rms error in hh:
  dhrms=sum((hh-wkp)**2)
   !Re-assign hh & htot:
  hh=wkp
  htot=one+hh

   !Compute relative vorticity zz = (1+hh)*(f+qq) - f:
  wkp=qq*htot
  qadd=-dsumi*sum(wkp)
   !qadd+qq is the corrected PV anomaly:
  qq=qq+qadd
  zz=htot*(cof+qq)-cof

   !Create spectral version of relative vorticity:
  call ptospc(ng,ng,zz,wkb,xfactors,yfactors,xtrig,ytrig)
   !wkb: vorticity in spectral space (zz destroyed but recreated below)

   !Solve Lap(wka) = zz spectrally:
  wka=rlap*wkb

   !Filter relative vorticity and bring back to physical space as zz:
  wkb=filt*wkb
  call spctop(ng,ng,wkb,zz,xfactors,yfactors,xtrig,ytrig)

   !Compute derivatives in spectral space:
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)
  wkb=-wkb

   !New velocity components in spectral space are in (wkb,wkd) here.

   !Convert to physical space as (hx,hy):
  call spctop(ng,ng,wkb,hx,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,hy,xfactors,yfactors,xtrig,ytrig)

   !Add mean flow (uio,vio):
  uio=-sum(hh*hx)*dsumi
  vio=-sum(hh*hy)*dsumi
  hx=hx+uio
  hy=hy+vio

   !Compute rms error in uu & vv:
  durms=sum((uu-hx)**2+(vv-hy)**2)

   !Re-assign velocity components:
  uu=hx
  vv=hy

   !Compute overall error:
  enorm=sqrt((durms+csq*dhrms)/sum(uu**2+vv**2+csq*hh**2+1.d-20))
enddo
 !Passing this, we have converged.

!-----------------------------------------------------------------
!Iterate to find the balanced fields:
hhpre=hh
uupre=uu
vvpre=vv

 !Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)
  !Obtain balanced estimate for gamma (gg):
  call jacob(uu,vv,wkp)
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  wkp=dd*uu
  wkq=dd*vv
  call divs(wkp,wkq,wka)
  gs=filt*(wka-two*wkb)
  wka=gs
  call spctop(ng,ng,wka,gg,xfactors,yfactors,xtrig,ytrig)

  !Obtain balanced estimate for delta (dd):
  wkp=hh*uu
  wkq=hh*vv
  call divs(wkp,wkq,wka)
  wkp=zz*uu
  wkq=zz*vv
  call divs(wkp,wkq,wkb)
  ds=helm*(cof*wkb-c2g2*wka)
  wka=ds
  call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)

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
   !Compute rms error in hh:
  dhrms=sum((hh-hhpre)**2)

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
   !Compute rms error in uu & vv:
  durms=sum((uu-uupre)**2+(vv-vvpre)**2)

   !Compute overall error:
  enorm=sqrt((durms+csq*dhrms)/sum(uupre**2+vvpre**2+csq*hhpre**2+1.d-20))

  write(*,*) ' Relative energy error = ',enorm

  !Otherwise continue with another iteration:
  hhpre=hh
  uupre=uu
  vvpre=vv
enddo

!-----------------------------------------------------------------
!Write data:

!PV anomaly (need to add Jacobian term):
call jacob(hh,dd,wkp)
call dealias(wkp)
qq=qq+hbsq3*wkp
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,qq
close(11)

!Divergence:
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,dd
close(11)

!Acceleration divergence:
open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,gg
close(11)

write(*,*)
write(*,*) ' Initial fields balanced and re-written.'

 !End main program
end program
!=======================================================================
