module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use constants
use sta2dfft

 !Common arrays, constants:
double precision:: rksq(ng,ng),opak(ng,ng),filt(ng,ng)
double precision:: c2g2(ng,ng),rlap(ng,ng),simp(ng,ng)
double precision:: bflo(ng,ng),bfhi(ng,ng),dissi(ng,ng)

 !For 2D FFTs:
double precision:: hrkx(ng),hrky(ng),rk(ng)
double precision:: xtrig(2*ng),ytrig(2*ng)
integer:: xfactors(5),yfactors(5)

double precision:: spmf(0:ng),alk(ng)
integer:: kmag(ng,ng),kmax,kmaxred


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================
subroutine init_spectral
! Initialises this module

implicit none

!Local variables:
double precision:: fac,rkmax,rks,snorm
double precision:: anu,rkfsq,fsq
integer:: kx,ky,k,kc

!----------------------------------------------------------------------
 !Set up 2D FFTs:
call init2dfft(ng,ng,twopi,twopi,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define wavenumbers and filtered wavenumbers:
rk(1)=zero
do k=1,ng/2-1
  rk(k+1)   =hrkx(2*k)
  rk(ng+1-k)=hrkx(2*k)
enddo
rk(ng/2+1)=hrkx(ng)

!-----------------------------------------------------------------------
 !Initialise arrays for computing the spectrum of any field:
rkmax=dble(ng/2)
kmax=nint(rkmax*sqrt(two))
do k=0,kmax
  spmf(k)=zero
enddo
do ky=1,ng
  do kx=1,ng
    k=nint(sqrt(rk(kx)**2+rk(ky)**2))
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(ng*ng) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi/dble(ng*ng)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=ng/2

!-----------------------------------------------------------------------
 !Define a variety of spectral operators:

 !Hyperviscosity coefficient (Dritschel, Gottwald & Oliver, JFM (2017)):
anu=cdamp*cof/rkmax**(2*nnu)
 !Assumes Burger number = 1.

 !Used for de-aliasing filter below:
rkfsq=(dble(ng)/3.d0)**2

 !Used for Butterworth filter below:
kc=ng/6
fac=one/dble(kc**2)

 !Squared Coriolis frequency:
fsq=cof**2

do ky=1,ng
  do kx=1,ng
    rks=rk(kx)**2+rk(ky)**2
     !Spectral c^2*grad^2 - f^2 operator:
    opak(kx,ky)=-(fsq+csq*rks)
     !Hyperviscous operator:
    dissi(kx,ky)=anu*rks**nnu
     !De-aliasing filter:
    if (rks .gt. rkfsq) then
      filt(kx,ky)=zero
      rksq(kx,ky)=zero
      c2g2(kx,ky)=zero
      rlap(kx,ky)=zero
      bflo(kx,ky)=zero
      bfhi(kx,ky)=zero
    else
      filt(kx,ky)=one
       !Squared wavenumber (-Laplace operator):
      rksq(kx,ky)=rks
       !c^2*grad^2:
      c2g2(kx,ky)=-csq*rks
       !grad^{-2}:
      rlap(kx,ky)=-one/(rks+small)
       !Butterworth low-pass (F) & high-pass (1-F) filters:
      bflo(kx,ky)=one/(one+(fac*rks)**2)
      bfhi(kx,ky)=one-bflo(kx,ky)
    endif
     !Semi-implicit operator for inverting divergence:
    simp(kx,ky)=one/((hdtisq-opak(kx,ky))*(one+dt2*dissi(kx,ky)))
     !Redfine damping operator for use in qd evolution:
    dissi(kx,ky)=one/(one+dt2*dissi(kx,ky))
  enddo
enddo

 !Ensure mean potentials and height anomaly remain zero:
rlap(1,1)=zero

return 
end subroutine

!======================================================================
subroutine main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
! Given the PV anomaly qs, divergence ds and acceleration divergence gs
! (all in spectral space), this routine computes the dimensionless depth 
! anomaly hh and the velocity field (uu,vv) in physical space.  It also 
! returns the corrected PV anomaly (qq) and relative vorticity (zz) in 
! physical space.

implicit none

 !Passed variables:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng) !Spectral
double precision:: hh(ng,ng),uu(ng,ng),vv(ng,ng) !Physical
double precision:: qq(ng,ng),zz(ng,ng)           !Physical

 !Local variables:
double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      hh, uu & vv from qq, dd & gg.  The energy error is computed from 
 !      <(u-u0)^2+(v-v0)^2+c^2*(h-h0)^2>/<u^2+v^2+c^2*h^2>
 !      where <:> means a domain average and (u0,v0,h0) is the previous
 !      guess for (u,v,h).

 !Physical work arrays:
double precision:: htot(ng,ng),hx(ng,ng),hy(ng,ng)
double precision:: gg(ng,ng),wkp(ng,ng)

 !Spectral work arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng)
double precision:: uds(ng,ng),vds(ng,ng)

 !Other constants:
double precision:: uio,vio
double precision:: qadd,qbar,fqbar
double precision:: dhrms,durms,enorm

!-------------------------------------------------------
 !Define total dimensionless height:
htot=one+hh

 !Define spectral divergent velocity (never changes in iteration):
wkc=rlap*ds
 !This solves Lap(wkc) = dd in spectral space
call xderiv(ng,ng,hrkx,wkc,uds)
call yderiv(ng,ng,hrky,wkc,vds)

 !Obtain a physical space copy of qs & gs:
wka=qs
call spctop(ng,ng,wka,qq,xfactors,yfactors,xtrig,ytrig)
wka=gs
call spctop(ng,ng,wka,gg,xfactors,yfactors,xtrig,ytrig)

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

   !Invert [c^2*grad^2-f*(f+<q>)]h = f*(q+h*(q-<q>) - gamma
  wkp=cof*(qq+hh*(qq-qbar))-gg
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

   !New velocity components in spectral space, written in (wkb,wkd):
  wkb=uds-wkb  !uds is the fixed divergent part of uu
  wkd=vds+wkd  !vds is the fixed divergent part of vv

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
  enorm=sqrt((durms+csq*dhrms)/sum(uu**2+vv**2+csq*hh**2))
enddo
 !Passing this, we have converged.

!------------------------------------------------------------------

return
end subroutine

!=================================================================
subroutine gradient(ff,ffx,ffy)
! Computes the gradient ffx = dF/dx & ffy = dF/dy of a field F.
! *** ff is in spectral space whereas (ffx,ffy) are in physical space

implicit none

 !Passed arrays:
double precision:: ff(ng,ng)             !Spectral
double precision:: ffx(ng,ng),ffy(ng,ng) !Physical

 !Local array:
double precision:: vtmp(ng,ng)           !Spectral

 !Get derivatives of F:
call xderiv(ng,ng,hrkx,ff,vtmp)
call spctop(ng,ng,vtmp,ffx,xfactors,yfactors,xtrig,ytrig)

call yderiv(ng,ng,hrky,ff,vtmp)
call spctop(ng,ng,vtmp,ffy,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!=================================================================
subroutine jacob(aa,bb,cc)
! Computes the Jacobian of aa and bb and returns it in cc.
! All passed variables are in physical space.

implicit none

 !Passed arrays:
double precision:: aa(ng,ng),bb(ng,ng),cc(ng,ng)           !Physical

 !Work arrays:
double precision:: ax(ng,ng),ay(ng,ng),bx(ng,ng),by(ng,ng) !Physical
double precision:: wka(ng,ng),wkb(ng,ng)                   !Spectral

!---------------------------------------------------------
cc=aa
call ptospc(ng,ng,cc,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
 !Get derivatives of aa:
call xderiv(ng,ng,hrkx,wka,wkb)
call spctop(ng,ng,wkb,ax,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,wkb)
call spctop(ng,ng,wkb,ay,xfactors,yfactors,xtrig,ytrig)

cc=bb
call ptospc(ng,ng,cc,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
 !Get derivatives of bb:
call xderiv(ng,ng,hrkx,wka,wkb)
call spctop(ng,ng,wkb,bx,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,wkb)
call spctop(ng,ng,wkb,by,xfactors,yfactors,xtrig,ytrig)

cc=ax*by-ay*bx

return
end subroutine

!=================================================================
subroutine divs(aa,bb,cs)
! Computes the divergence of (aa,bb) and returns it in cs.
! Both aa and bb in physical space but cs is in spectral space.

implicit none

 !Passed arrays:
double precision:: aa(ng,ng),bb(ng,ng)   !Physical
double precision:: cs(ng,ng)             !Spectral

 !Work arrays:
double precision:: wkp(ng,ng)            !Physical
double precision:: wka(ng,ng),wkb(ng,ng) !Spectral

!---------------------------------------------------------
wkp=aa
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call xderiv(ng,ng,hrkx,wka,wkb)

wkp=bb
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,cs)

cs=wkb+cs

return
end subroutine

!===================================================================

subroutine spec1d(ss,spec)
! Computes the 1d spectrum of a spectral field ss and returns the
! result in spec.

implicit none

 !Passed variables:
double precision:: ss(ng,ng),spec(0:ng)

 !Local variables:
integer:: kx,ky,k

!--------------------------------------------------------
do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(1,1)
spec(k)=spec(k)+f14*ss(1,1)**2

 !y-independent mode:
do kx=2,ng
  k=kmag(kx,1)
  spec(k)=spec(k)+f12*ss(kx,1)**2
enddo

 !x-independent mode:
do ky=2,ng
  k=kmag(1,ky)
  spec(k)=spec(k)+f12*ss(1,ky)**2
enddo

 !All other modes:
do ky=2,ng
  do kx=2,ng
    k=kmag(kx,ky)
    spec(k)=spec(k)+ss(kx,ky)**2
  enddo
enddo

return
end subroutine

!===================================================================

end module     
