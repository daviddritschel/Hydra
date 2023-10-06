module spectral

use constants
use sta2dfft

 !Maximum x & y wavenumbers:
integer,parameter:: nwx=nx/2,nwy=ny/2

 !Common arrays, constants:
double precision:: rksq(nx,ny),qdiss(nx,ny),adiss(nx,ny)
double precision:: green(nx,ny),filt(nx,ny),flo(nx,ny),fhi(nx,ny)
double precision:: hrkx(nx),hrky(ny)
double precision:: bety(ny)

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(nx,ny),kmax,kmaxred

!=======================================================================!
! From main code: call init_spectral             to initialise          !
! then            call main_invert(qq,uu,vv,pp)  to perform inversion   !
!=======================================================================!

contains

!=================================================================

subroutine init_spectral

implicit none

 !Local variables:
double precision:: rkx(nx),rky(ny)
double precision:: scx,rkxmax,scy,rkymax,rkfsq,rkmsi
double precision:: delk,delki,fac,rksqu,snorm
integer, dimension(:), allocatable :: seed
integer:: kx,ky,k,iy,i

!----------------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define x wavenumbers:
rkx(1)=zero
do kx=1,nwx-1
  rkx(kx+1)   =hrkx(2*kx)
  rkx(nx+1-kx)=hrkx(2*kx)
enddo
rkx(nwx+1)=hrkx(nx)

scx=twopi/ellx
rkxmax=scx*dble(nwx)

 !Define y wavenumbers:
rky(1)=zero
do ky=1,nwy-1
  rky(ky+1)   =hrky(2*ky)
  rky(ny+1-ky)=hrky(2*ky)
enddo
rky(nwy+1)=hrky(ny)

scy=twopi/elly
rkymax=scy*dble(nwy)

!-----------------------------------------------------------------------
 !Define Green function:
green(1,1)=zero
do kx=2,nx
  green(kx,1)=-one/(rkx(kx)**2+kdsq)
enddo
do ky=2,ny
  green(1,ky)=-one/(rky(ky)**2+kdsq)
enddo
do ky=2,ny
  do kx=2,nx
    green(kx,ky)=-one/(rkx(kx)**2+rky(ky)**2+kdsq)
  enddo
enddo

!-----------------------------------------------------------------------
 !Define squared total wavenumber, diffusion operators, low and hi-pass
 !Butterworth filters (Wireless Engineer, vol. 7, 1930, pp. 536-541),
 !hyperdiffusion for qd, and the de-aliasing filter:
fac=36.d0/dble(ngridp)
rkfsq=dble(ngridp)/9.d0
rkmsi=one/max(rkxmax**2,rkymax**2)
do ky=1,ny
  do kx=1,nx
    rksqu=rkx(kx)**2+rky(ky)**2
    qdiss(kx,ky)=cdamp*(rkmsi*rksqu)**nnu
    adiss(kx,ky)=eta*rksqu
    if (rksqu .gt. rkfsq) then
      filt(kx,ky)=zero
      flo(kx,ky)=zero
      fhi(kx,ky)=zero
      rksq(kx,ky)=zero
      green(kx,ky)=zero
    else
      filt(kx,ky)=one
      flo(kx,ky)=one/(one+(fac*rksqu)**2)
      fhi(kx,ky)=one-flo(kx,ky)
      rksq(kx,ky)=rksqu
    endif
  enddo
enddo

!-----------------------------------------------------------------------
 !Initialise arrays for computing the spectrum of any field:
delk=sqrt(scx*scy)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do ky=1,ny
  do kx=1,nx
    k=nint(sqrt(rkx(kx)**2+rky(ky)**2)*delki)
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(nx*ny) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi*dsumi   !dsumi = 1/dble(nx*ny)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(delk*dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=nint(sqrt((rkxmax**2+rkymax**2)/two)*delki)

if (beffect) then
   !Define beta*y if beta is non-zero:
  do iy=1,ny
    bety(iy)=beta*(gly*dble(iy-1)+ymin)
  enddo
endif

!-----------------------------------------------------------------------
 !Initialise random number generator if forcing is present:
if (zforcing .or. aforcing) then
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:)=iseed
  do i=1,iseed
    call random_seed(put=seed)
  enddo
endif

return 
end subroutine init_spectral

!=================================================================

subroutine main_invert(qq,uu,vv,pp)
! Given the PV anomaly qq in spectral space, this routine computes
! the streamfunction pp in spectral space and the velocity field
! (uu,vv) in physical space.

implicit none

 !Passed variables:
double precision:: qq(nx,ny),pp(nx,ny) !Spectral
double precision:: uu(ny,nx),vv(ny,nx) !Physical

 !Local variables:
double precision:: vtmp(nx,ny)

!--------------------------------------------------------
 !Solve for psi (pp) by inverting Lap - k_d^2 spectrally:
pp=green*qq

 !Get velocity field:
call xderiv(nx,ny,hrkx,pp,vtmp)
call spctop(nx,ny,vtmp,vv,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,pp,vtmp)
call spctop(nx,ny,vtmp,uu,xfactors,yfactors,xtrig,ytrig)

 !Copy -uu into uu:
uu=-uu
 !(green is spectrally truncated => uu & vv are as well)

return
end subroutine main_invert

!=================================================================

subroutine lorentz(apot,aax,aay,clf)
! Computes the magnetic field induced part of the PV source (clf),
!      S_b = B_0*dj/dx - J(A,j)
! where j = -Lap{A} is the current density.
! Note: apot = A in *** spectral space *** but all other fields
!       are in physical space.

implicit none

 !Passed arrays:
double precision:: apot(nx,ny) !Spectral (input)
double precision:: aax(ny,nx),aay(ny,nx),clf(ny,nx) !Physical (output)

 !Local variables:
double precision:: jjx(ny,nx),jjy(ny,nx) !Physical
double precision:: jj(nx,ny),vtmp(nx,ny) !Spectral

!---------------------------------------------------
 !Get derivatives of A:
call xderiv(nx,ny,hrkx,apot,vtmp)
call spctop(nx,ny,vtmp,aax,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,apot,vtmp)
call spctop(nx,ny,vtmp,aay,xfactors,yfactors,xtrig,ytrig)

 !Define j in spectral space:
jj=rksq*apot

 !Get derivatives of j:
call xderiv(nx,ny,hrkx,jj,vtmp)
call spctop(nx,ny,vtmp,jjx,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,jj,vtmp)
call spctop(nx,ny,vtmp,jjy,xfactors,yfactors,xtrig,ytrig)

 !Compute source term, S_b:
clf=(b0+aay)*jjx-aax*jjy

return
end subroutine lorentz

!===============================================================

subroutine gradient(ff,ffx,ffy)
! Computes the gradient ffx = dF/dx & ffy = dF/dy of a field F.
! *** ff is in spectral space whereas (ffx,ffy) are in physical space

implicit none

 !Passed arrays:
double precision:: ff(nx,ny) !Spectral (note order)
double precision:: ffx(ny,nx),ffy(ny,nx) !Physical

 !Local array:
double precision:: vtmp(nx,ny)

 !Get derivatives of F:
call xderiv(nx,ny,hrkx,ff,vtmp)
call spctop(nx,ny,vtmp,ffx,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,ff,vtmp)
call spctop(nx,ny,vtmp,ffy,xfactors,yfactors,xtrig,ytrig)

return
end subroutine gradient

!===================================================================

subroutine spec1d(ss,spec,iopt)
! Computes the 1d spectrum of a spectral field ss and returns the
! result in spec.
! If iopt = 1, the spectral field is multiplied first by K^2
!              and is returned modified!

implicit none

 !Passed variables:
double precision:: ss(nx,ny),spec(0:max(nx,ny))
integer:: iopt

 !Local variables:
integer:: kx,ky,k

!--------------------------------------------------------
if (iopt .eq. 1) then
   !Multiply spectral field by K^2 (filtered, see init_spectral):
  ss=ss*rksq
endif

do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(1,1)
spec(k)=spec(k)+f14*ss(1,1)**2

 !y-independent mode:
do kx=2,nx
  k=kmag(kx,1)
  spec(k)=spec(k)+f12*ss(kx,1)**2
enddo

 !x-independent mode:
do ky=2,ny
  k=kmag(1,ky)
  spec(k)=spec(k)+f12*ss(1,ky)**2
enddo

 !All other modes:
do ky=2,ny
  do kx=2,nx
    k=kmag(kx,ky)
    spec(k)=spec(k)+ss(kx,ky)**2
  enddo
enddo

return
end subroutine spec1d

!===================================================================

subroutine ranspec(ff,rms,pow,k0)
! Computes a random field ff with spectrum c k^{2p-1} exp[-2(k/k_0)^2]
! where p = pow, k_0 = k0 and c is determined by the rms value of the
! resulting physical field.  ff is returned in physical space.

implicit none

 !Passed variables:
double precision:: ff(ny,nx),rms,pow
integer:: k0

 !Local variables:
double precision:: ss(nx,ny)
double precision:: fac,p1,s,uni,phix,phiy
double precision:: amp,cx,sx,cy,sy
integer:: kx,ky,kxc,kyc

!--------------------------------------------------------
! Generate spectrum / k (actually, its square root):
fac=one/dble(k0*k0)
if (pow .gt. 1.000001d0) then
  p1=pow-one
  do ky=1,nwy+1
    do kx=1,nwx+1
      s=fac*rksq(kx,ky)
      ss(kx,ky)=sqrt(fac*s**p1*exp(-two*s))
    enddo
  enddo
else
   !Here pow = 1 is assumed:
  do ky=1,nwy+1
    do kx=1,nwx+1
      ss(kx,ky)=filt(kx,ky)*exp(-fac*rksq(kx,ky))
    enddo
  enddo
endif

! Apply to generate full spectrum:
do ky=2,nwy
  kyc=ny+2-ky
  do kx=2,nwx
    kxc=nx+2-kx
    call random_number(uni)
    phix=twopi*uni-pi
    call random_number(uni)
    phiy=twopi*uni-pi
    cx=cos(phix)
    sx=sin(phix)
    cy=cos(phiy)
    sy=sin(phiy)
    amp=ss(kx,ky)
    ss(kx ,ky )=amp*cx*cy
    ss(kxc,ky )=amp*sx*cy
    ss(kx, kyc)=amp*cx*sy
    ss(kxc,kyc)=amp*sx*sy
  enddo
enddo

ky=1
do kx=2,nwx
  kxc=nx+2-kx
  call random_number(uni)
  phix=twopi*uni-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cx
  ss(kxc,ky )=amp*sx
enddo

kx=1
do ky=2,nwy
  kyc=ny+2-ky
  call random_number(uni)
  phiy=twopi*uni-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cy
  ss(kx, kyc)=amp*sy
enddo

ky=nwy+1
do kx=2,nwx
  kxc=nx+2-kx
  call random_number(uni)
  phix=twopi*uni-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cx
  ss(kxc,ky )=amp*sx
enddo

kx=nwx+1
do ky=2,nwy
  kyc=ny+2-ky
  call random_number(uni)
  phiy=twopi*uni-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cy
  ss(kx, kyc)=amp*sy
enddo

ss(1,1)=zero
ss(nwx+1,nwy+1)=zero

! Transform to physical space:
call spctop(nx,ny,ss,ff,xfactors,yfactors,xtrig,ytrig)

! Work out rms:
fac=sqrt(dsumi*sum(ff**2))   !dsumi = 1/dble(nx*ny)

! Renormalise field:
fac=rms/fac
ff=fac*ff

return
end subroutine ranspec

!=======================================================

end module
