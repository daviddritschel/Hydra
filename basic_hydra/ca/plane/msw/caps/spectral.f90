module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use constants
use sta2dfft

 !Common arrays, constants:
double precision:: rksq(ng,ng),rlap(ng,ng)
double precision:: opak(ng,ng),gwop(ng,ng),helm(ng,ng)
double precision:: bflo(ng,ng),bfhi(ng,ng),filt(ng,ng)
double precision:: qdis(ng,ng),bdis(ng,ng)

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

 !Hyperviscosity coefficient (except for zeta_rms multiplier):
anu=cdamp/rkmax**(2*nnu)

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
     !Hyperviscous operator (except for zeta_rms multiplier):
    qdis(kx,ky)=anu*rks**nnu
     !Magnetic diffusivity operator:
    bdis(kx,ky)=eta*rks
     !De-aliasing filter:
    if (rks .gt. rkfsq) then
      filt(kx,ky)=zero
      rksq(kx,ky)=zero
      gwop(kx,ky)=zero
      helm(kx,ky)=zero
      rlap(kx,ky)=zero
      bflo(kx,ky)=zero
      bfhi(kx,ky)=zero
    else
      filt(kx,ky)=one
       !Squared wavenumber (-Laplace operator):
      rksq(kx,ky)=rks
       !G = c^2*grad^2 - f^2 (filtered):
      gwop(kx,ky)=opak(kx,ky)
       !H = G^{-1} (filtered):
      helm(kx,ky)=one/opak(kx,ky)
       !grad^{-2}:
      rlap(kx,ky)=-one/(rks+small)
       !Butterworth low-pass (F) & high-pass (1-F) filters:
      bflo(kx,ky)=one/(one+(fac*rks)**2)
      bfhi(kx,ky)=one-bflo(kx,ky)
    endif
  enddo
enddo

 !Ensure mean potentials and height anomaly remain zero:
rlap(1,1)=zero

return 
end subroutine

!======================================================================
subroutine main_invert(qs,ds,gs,hh,uu,vv,zz)
! Given the linearised PV anomaly qs, divergence ds and acceleration
! divergence gs (all in spectral space), this routine computes the
! dimensionless height anomaly hh, the velocity (uu,vv) and the
! relative (vertical) vorticity zz in physical space.

implicit none

 !Passed variables:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)           !Spectral
double precision:: hh(ng,ng),uu(ng,ng),vv(ng,ng),zz(ng,ng) !Physical

 !Local variables:

 !Spectral work arrays:
double precision:: wka(ng,ng),wkb(ng,ng)
double precision:: wkc(ng,ng),wkd(ng,ng),wke(ng,ng)

 !Other constants:
double precision:: uio,vio

!------------------------------------------------------------
 !Obtain dimensionless height anomaly and relative vorticity:
wka=helm*(cof*qs-gs) ! spectral h_tilde
wkc=qs+cof*wka       ! spectral zeta
call spctop(ng,ng,wka,hh,xfactors,yfactors,xtrig,ytrig)
 !hh now contains h_tilde on the grid

 !Define spectral non-divergent velocity:
wke=rlap*wkc
 !This solves Lap(wke) = zeta in spectral space
call xderiv(ng,ng,hrkx,wke,wka) ! spectral dpsi/dx
call yderiv(ng,ng,hrky,wke,wkb) ! spectral dpsi/dy

call spctop(ng,ng,wkc,zz,xfactors,yfactors,xtrig,ytrig)
 !zz now contains zeta on the grid

 !Define spectral divergent velocity:
wke=rlap*ds
 !This solves Lap(wke) = dd in spectral space
call xderiv(ng,ng,hrkx,wke,wkc) ! spectral dxi/dx
call yderiv(ng,ng,hrky,wke,wkd) ! spectral dxi/dy

 !Add up components to define velocity field:
wkc=wkc-wkb
wkd=wkd+wka
call spctop(ng,ng,wkc,uu,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,wkd,vv,xfactors,yfactors,xtrig,ytrig)

 !Add mean flow (uio,vio) consistent with zero momentum:
uio=-sum(hh*uu)*dsumi
vio=-sum(hh*vv)*dsumi
uu=uu+uio
vv=vv+vio

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
