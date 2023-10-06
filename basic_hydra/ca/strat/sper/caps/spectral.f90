module spectral

use constants
use sta2dfft

 !Declarations:
implicit none

 !Common arrays, constants:
double precision:: yh0(0:ny),yh1(0:ny),pbar(0:ny)
double precision:: rkx(0:nxm1),hrkx(nx),rky(ny)
double precision:: green(0:nxm1,0:ny),hdis(0:nxm1,0:ny)
double precision::  bflo(0:nxm1,0:ny),bfhi(0:nxm1,0:ny)
double precision::  filt(0:nxm1,0:ny),decy(nym1,0:nxm1)
double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)
integer:: kmag(0:nxm1,0:ny),kmax

!============================================================================!
! From main code: call init_spectral                  to initialise          !
! then            call main_invert(zs,uu,vv,zz,zavg)  to perform inversion   !
!============================================================================!

contains

!=====================================================================

subroutine init_spectral

 !Declarations:
implicit none

 !Local variables:

 !Maximum x wavenumber:
integer,parameter:: nwx=nx/2

 !Others:
double precision:: fac,yg,scx,scy,rkxmax,rkymax,rkfsq
double precision:: delk,delki,snorm,div,visc
integer:: iy,kx,kxc,ky,k

!---------------------------------------------------------------------
 !Fractional y grid values: 
fac=one/dble(ny)
do iy=0,ny
  yh1(iy)=fac*dble(iy)
  yh0(iy)=one-yh1(iy)
enddo

 !Define part of streamfunction proportional to the mean vorticity:
do iy=0,ny
  yg=gly*dble(iy)
  pbar(iy)=f12*yg*(yg-elly)
enddo

!---------------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,rky)

 !Define x wavenumbers:
rkx(0)=zero
do kx=1,nwx-1
  kxc=nx-kx
  rkx(kx )=hrkx(2*kx)
  rkx(kxc)=hrkx(2*kx)
enddo
rkx(nwx)=hrkx(nx)

 !Initialise arrays for computing the spectrum of any field:
scx=twopi/ellx
rkxmax=scx*dble(nwx)
scy=pi/elly
rkymax=scy*dble(ny)
delk=sqrt(scx**2+scy**2)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do kx=0,nxm1
  k=nint(rkx(kx)*delki)
  kmag(kx,ky)=k
  spmf(k)=spmf(k)+one
enddo
do ky=1,ny
  do kx=0,nxm1
    k=nint(sqrt(rkx(kx)**2+rky(ky)**2)*delki)
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(nx*ny) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi*dsumi !dsumi = 1/dble(nx*ny)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(delk*dble(k))
enddo

!---------------------------------------------------------------------
 !Define de-aliasing filter (2/3 rule):
rkfsq=(f23*rkxmax)**2
filt(0,0)=one
do kx=1,nxm1
  if (rkx(kx)**2 .gt. rkfsq) then
    filt(kx,0)=zero
  else
    filt(kx,0)=one
  endif
enddo
do ky=1,ny
  if (rky(ky)**2 .gt. rkfsq) then
    filt(0,ky)=zero
  else
    filt(0,ky)=one
  endif
enddo
do ky=1,ny
  do kx=1,nxm1
    if (rkx(kx)**2+rky(ky)**2 .gt. rkfsq) then
      filt(kx,ky)=zero
    else
      filt(kx,ky)=one
    endif
  enddo
enddo

!---------------------------------------------------------------------
 !Define Butterworth filter:
fac=9.d0/rkxmax**2
bflo(0,0)=one
do kx=1,nxm1
  bflo(kx,0)=one/(one+(fac*rkx(kx)**2)**2)
enddo
do ky=1,ny
  bflo(0,ky)=one/(one+(fac*rky(ky)**2)**2)
enddo
do ky=1,ny
  do kx=1,nxm1
    bflo(kx,ky)=one/(one+(fac*(rkx(kx)**2+rky(ky)**2))**2)
  enddo
enddo
bfhi=filt*(one-bflo)
bflo=filt*bflo

!---------------------------------------------------------------------
 !Define Green function:
green(0,0)=zero
do kx=1,nxm1
  green(kx,0)=-one/rkx(kx)**2
enddo
do ky=1,ny
  green(0,ky)=-one/rky(ky)**2
enddo
do ky=1,ny
  do kx=1,nxm1
    green(kx,ky)=-one/(rkx(kx)**2+rky(ky)**2)
  enddo
enddo

!---------------------------------------------------------------------
 !Hyperbolic functions used for solutions of Laplace's equation:
decy(1:nym1,0)=yh1(1:nym1)
do kx=1,nxm1
  fac=rkx(kx)*elly
  div=one/(one-exp(-two*fac))
  decy(1:nym1,kx)=(exp(-fac*(one-yh1(1:nym1)))- &
                   exp(-fac*(one+yh1(1:nym1))))*div
enddo

!---------------------------------------------------------------------
 !Define dissipation operator:
visc=prediss/max(rkxmax,rkymax)**(2*nnu)
hdis(0,0)=zero
do kx=1,nxm1
  hdis(kx,0)=visc*rkx(kx)**(2*nnu)
enddo
do ky=1,ny
  hdis(0,ky)=visc*rky(ky)**(2*nnu)
enddo
do ky=1,ny
  do kx=1,nxm1
    hdis(kx,ky)=visc*(rkx(kx)**2+rky(ky)**2)**nnu
  enddo
enddo

return 
end subroutine init_spectral

!=====================================================================

subroutine main_invert(zs,uu,vv,zz,zavg)
! Given the vorticity zs in spectral space, and its average value zavg
! in physical space, this routine computes the velocity (uu,vv) and
! the vorticity zz in physical space.
  
 !Declarations:
implicit none

 !Input array (spectral):
double precision:: zs(0:nxm1,0:ny)
 !Input average vorticity:
double precision:: zavg

!Output arrays (physical):
double precision:: zz(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)

 !Local arrays:
double precision:: ss(0:nxm1,0:ny),pp(0:ny,0:nxm1)
double precision:: pbot(0:nxm1),ptop(0:nxm1),cppy(nym1,0:nxm1)

 !Other quantities:
integer:: kx,ky,ix,iy

!--------------------------------
 !Solve for psi (pp):

 !(1) Compute vorticity zz in physical space:
ss=zs !***warning*** zs(0,0) = 0 is assumed
call spctop_fc(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)
zz=zz+zavg
 
 !(2) Invert spectral vorticity zs to get uncorrected streamfunction pp:
ss=green*zs
call spctop_fc(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

 !(3) Add part of pp due to mean vorticity zavg:
do ix=0,nxm1
  pp(:,ix)=pp(:,ix)+zavg*pbar
enddo

 !(4) Do a sine transform of pp at y = ymin and ymax and obtain the
 !    interior field (cppy) that must be subtracted to give pp = 0
 !    at y = ymin and ymax:
pbot=pp(0,:)
ptop=pp(ny,:)
call forfft(1,nx,pbot,xtrig,xfactors)
call forfft(1,nx,ptop,xtrig,xfactors)

 !Define the interior semi-spectral field:
do kx=0,nxm1
  do iy=1,nym1
    cppy(iy,kx)=pbot(kx)*decy(ny-iy,kx)+ptop(kx)*decy(iy,kx)
  enddo
enddo
 !Invert using a full transform in x:
call revfft(nym1,nx,cppy,xtrig,xfactors)

 !(5) Remove cppy to obtain the final streamfunction pp:
pp(0,:)=zero
pp(ny,:)=zero
pp(1:nym1,:)=pp(1:nym1,:)-cppy(1:nym1,:)

 !(6) Compute velocity field from pp:
call getvel(pp,uu,vv)

return
end subroutine main_invert

!=====================================================================

subroutine getvel(pp,uu,vv)
! Computes the velocity components uu & vv from the streamfunction
! pp via uu = -d(pp)/dy and vv = d(pp)/dx.
! *** pp, uu & vv are all in physical space
! *** and include the domain edges.

implicit none

 !Passed arrays:
double precision:: pp(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)

 !Local arrays:
double precision:: ppi(ny  ,0:nxm1),pps(0:nxm1  ,ny)
double precision:: ppx(0:nxm1,  ny),vvi(  ny,0:nxm1)
double precision:: ppy(0:nxm1,0:ny)

!-------------------------------------------------------------------
 !Copy non-zero interior values of pp to ppi:
ppi(1:nym1,:)=pp(1:nym1,:)

 !Transform ppi to spectral space:
call ptospc_fs(nx,ny,ppi,pps,xfactors,yfactors,xtrig,ytrig)

 !Apply de-aliasing filter:
pps(:,1:ny)=pps(:,1:ny)*filt(:,1:ny)

 !Compute d(ppi)/dx = ppx spectrally:
call xderiv_fs(nx,ny,hrkx,pps,ppx)

 !Transform ppx back to physical space as vvi:
call spctop_fs(nx,ny,ppx,vvi,xfactors,yfactors,xtrig,ytrig)

 !Copy vvi into vv and add on zero edge values at iy = 0 & ny:
vv(0,:)=zero
vv(1:nym1,:)=vvi(1:nym1,:)
vv(ny,:)=zero

!-------------------------------------------------------------------
 !Compute d(ppi)/dy = ppy spectrally:
call yderiv_fs(nx,ny,rky,pps,ppy)

 !Transform ppy back to physical space as uu:
call spctop_fc(nx,ny,ppy,uu,xfactors,yfactors,xtrig,ytrig)

 !Correct sign:
uu=-uu

return
end subroutine getvel

!=====================================================================

subroutine spec1d_fc(var,spec)
! Computes the 1d spectrum of a spectral field var which is
! periodic in x and represented by a cosine series in y.
! Returns the result in spec.

implicit none

 !Passed arrays:
double precision:: var(0:nxm1,0:ny),spec(0:max(nx,ny))
 !Local indices:
integer:: k,kx,ky

!-------------------------------------------------------------------
 !Initialise spectrum:
spec(0:kmax)=zero

 !x and y-independent mode:
k=kmag(0,0)
spec(k)=spec(k)+f14*var(0,0)**2

 !y-independent mode:
do kx=1,nxm1
  k=kmag(kx,0)
  spec(k)=spec(k)+f12*var(kx,0)**2
enddo

 !x-independent mode:
do ky=1,ny
  k=kmag(0,ky)
  spec(k)=spec(k)+f12*var(0,ky)**2
enddo

 !All other modes:
do ky=1,ny
  do kx=1,nxm1
    k=kmag(kx,ky)
    spec(k)=spec(k)+var(kx,ky)**2
  enddo
enddo

return
end subroutine spec1d_fc

!=======================================================================

end module spectral
