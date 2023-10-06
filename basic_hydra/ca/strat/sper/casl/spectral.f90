module spectral

use constants
use sta2dfft

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Maximum x wavenumber:
integer,parameter:: nwx=nx/2,nwxm1=nwx-1,nwxp1=nwx+1

 !Common arrays, constants:
double precision:: pbar(0:ny)
double precision:: decy(nym1,0:nxm1)
double precision:: yh0(0:ny),yh1(0:ny)
double precision:: green(0:nxm1,0:ny)
double precision:: rkx(0:nxm1),hrkx(nx),rky(ny)

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

double precision:: frkx(0:nxm1),frky(ny)
double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(0:nxm1,0:ny),kmax,kmaxred,ifail

!====================================================================!
! From main code: call init_invert            to initialise          !
! then            call main_invert(zz,uu,vv)  to perform inversion   !
!====================================================================!

contains

!===========================
subroutine init_spectral

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

!----------------------

 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,rky)

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

 !Define x wavenumbers:
rkx(0)=zero
do kx=1,nwxm1
  kxc=nx-kx
  rkx(kx )=hrkx(2*kx)
  rkx(kxc)=hrkx(2*kx)
enddo
rkx(nwx)=hrkx(nx)

scx=twopi/ellx
rkxmax=scx*dble(nwx)
frkx(0)=zero
do kx=1,nxm1
  wratx=rkx(kx)/rkxmax
  frkx(kx)=rkx(kx)*exp(-36.d0*wratx**36.d0)
enddo

 !Define y wavenumbers:
scy=pi/elly
rkymax=scy*dble(ny)
do ky=1,ny
  wraty=rky(ky)/rkymax
  frky(ky)=rky(ky)*exp(-36.d0*wraty**36.d0)
enddo
 
 !Initialise arrays for computing the spectrum of any field:
delk=sqrt(scx**2+scy**2)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do ky=0,ny
  do kx=0,nxm1
    k=nint(sqrt((scx*dble(kx))**2+(scy*dble(ky))**2)*delki)
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(nx*ny) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi/dble(nx*ny)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(delk*dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=nint(sqrt((rkxmax**2+rkymax**2)/two)*delki)

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

 !Hyperbolic functions used for solutions of Laplace's equation:
do iy=1,nym1
  decy(iy,0)=yh1(iy)
enddo
do kx=1,nxm1
  fac=rkx(kx)*elly
  div=one/(one-exp(-two*fac))
  do iy=1,nym1
    argm=fac*(one-yh1(iy))
    argp=fac*(one+yh1(iy))
    decy(iy,kx)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

return 
end subroutine

!===========================
subroutine main_invert(zz,uu,vv)

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: zz(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
 !Local arrays:
double precision:: ss(0:nxm1,0:ny),pp(0:ny,0:nxm1)
double precision:: pbot(0:nxm1),ptop(0:nxm1),cppy(nym1,0:nxm1)

!--------------------------------
 !Solve for psi (pp):

 !(1) compute mean zz (zbar):
zbar=zero
do ix=0,nxm1
  zbar=zbar+zz(0,ix)+zz(ny,ix)
enddo
zbar=f12*zbar
do ix=0,nxm1
  do iy=1,nym1
    zbar=zbar+zz(iy,ix)
  enddo
enddo
zbar=zbar/dble(nx*ny)

 !(2) Remove mean vorticity from zz:
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix)=zz(iy,ix)-zbar
  enddo
enddo
 
 !(3) FFT zz and invert to get uncorrected streamfunction pp:
call ptospc_fc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)
do ky=0,ny
  do kx=0,nxm1
    ss(kx,ky)=green(kx,ky)*ss(kx,ky)
  enddo
enddo
call spctop_fc(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

 !(4) Add part of pp due to mean vorticity:
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix)=pp(iy,ix)+zbar*pbar(iy)
  enddo
enddo

 !(5) Do a sine transform of pp at y = ymin and ymax and obtain the
 !    interior field (cppy) that must be subtracted to give pp = 0
 !    at y = ymin and ymax:
do ix=0,nxm1
  pbot(ix)=pp(0,ix)
  ptop(ix)=pp(ny,ix)
enddo
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

 !(6) Remove cppy to obtain the final streamfunction pp:
do ix=0,nxm1
  pp(0, ix)=zero
  pp(ny,ix)=zero
enddo

do ix=0,nxm1
  do iy=1,nym1
    pp(iy,ix)=pp(iy,ix)-cppy(iy,ix)
  enddo
enddo

call getvel(pp,uu,vv)

return
end subroutine

!===================================================================

subroutine getvel(pp,uu,vv)
! Computes the velocity components uu & vv from the streamfunction
! pp via uu = -d(pp)/dy and vv = d(pp)/dx.
! *** pp, uu & vv are all in physical space
! *** and include the domain edges.

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: pp(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)

 !Local arrays:
double precision:: ppi(ny  ,0:nxm1),pps(0:nxm1  ,ny)
double precision:: ppx(0:nxm1,  ny),vvi(  ny,0:nxm1)
double precision:: ppy(0:nxm1,0:ny),uui(0:ny,0:nxm1)

 !Copy non-zero interior values of pp to ppi:
do ix=0,nxm1
  do iy=1,nym1
    ppi(iy,ix)=pp(iy,ix)
  enddo
enddo

 !Transform ppi to spectral space:
call ptospc_fs(nx,ny,ppi,pps,xfactors,yfactors,xtrig,ytrig)

 !Compute d(ppi)/dx = ppx spectrally:
call xderiv_fs(nx,ny,hrkx,pps,ppx)

 !Transform ppx back to physical space as vvi:
call spctop_fs(nx,ny,ppx,vvi,xfactors,yfactors,xtrig,ytrig)

 !Copy vvi into vv and add on zero edge values at iy = 0 & ny:
do ix=0,nxm1
  vv(0,ix)=zero
  do iy=1,nym1
    vv(iy,ix)=vvi(iy,ix)
  enddo
  vv(ny,ix)=zero
enddo

 !Compute d(ppi)/dy = ppy spectrally:
call yderiv_fs(nx,ny,rky,pps,ppy)

 !Transform ppy back to physical space as uui:
call spctop_fc(nx,ny,ppy,uui,xfactors,yfactors,xtrig,ytrig)

 !Copy -uui into uu and add mean reference velocity:
do ix=0,nxm1
  do iy=0,ny
    uu(iy,ix)=uref-uui(iy,ix)
  enddo
enddo

return
end subroutine

!===================================================================

subroutine spec1d_fc(rvar,spec)
! Computes the 1d spectrum of a field rvar which is
! periodic in x and represented by a cosine series in y
! and returns the result in spec.
!   *** Warning: rvar is modified by this routine.***

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(0:ny,0:nxm1),spec(0:max(nx,ny))
 !Local array:
double precision:: wks(0:nxm1,0:ny)

 !Transform rvar to spectral space:
call ptospc_fc(nx,ny,rvar,wks,xfactors,yfactors,xtrig,ytrig)

do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(0,0)
spec(k)=spec(k)+f14*wks(0,0)**2

 !y-independent mode:
do kx=1,nxm1
  k=kmag(kx,0)
  spec(k)=spec(k)+f12*wks(kx,0)**2
enddo

 !x-independent mode:
do ky=1,ny
  k=kmag(0,ky)
  spec(k)=spec(k)+f12*wks(0,ky)**2
enddo

 !All other modes:
do ky=1,ny
  do kx=1,nxm1
    k=kmag(kx,ky)
    spec(k)=spec(k)+wks(kx,ky)**2
  enddo
enddo

return
end subroutine

!===================================================================
subroutine bxderiv(bb,bbx)

! Subroutine to take the x derivative of a field bb, assumed to
! be even across the x boundaries.


implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: bb(0:ny,0:nxm1),bbx(0:ny,0:nxm1)

 !Copy bb into bbx for FFT:
do ix=0,nxm1
  do iy=0,ny
    bbx(iy,ix)=bb(iy,ix)
  enddo
enddo

 !Carry out a full x transform on bbx:
call forfft(nyp1,nx,bbx,xtrig,xfactors)

 !Take derivative spectrally by wavenumber multiplication:
do iy=0,ny
  bbx(iy,  0)=zero
  bbx(iy,nwx)=zero
enddo

do kx=1,nwxm1
  fac=frkx(kx)
  kxc=nx-kx
  do iy=0,ny
    bxtmp      =-fac*bbx(iy,kxc)
    bbx(iy,kxc)= fac*bbx(iy, kx)
    bbx(iy, kx)= bxtmp
  enddo
enddo

 !Carry out a full inverse x transform on bbx:
call revfft(nyp1,nx,bbx,xtrig,xfactors)

return
end subroutine

!============================================================
subroutine byderiv(bb,bby)

! Subroutine to take the y derivative of a field bb, assumed to
! be even across the y boundaries.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: bb(0:ny,0:nxm1),bby(0:ny,0:nxm1)

 !Local array:
double precision:: b1d(0:ny),by1d(ny)

 !Loop over columns, computing db/dy along each x grid line:
do ix=0,nxm1
   !Copy bb into b1d for FFT:
  do iy=0,ny
    b1d(iy)=bb(iy,ix)
  enddo

   !Carry out a y cosine transform on b1d:
  call dct(1,ny,b1d,ytrig,yfactors)

   !Take derivative spectrally by wavenumber multiplication:
  do ky=1,nym1
    by1d(ky)=-frky(ky)*b1d(ky)
  enddo
  by1d(ny)=zero

   !Carry out a y sine transform on by1d:
  call dst(1,ny,by1d,ytrig,yfactors)

   !Copy by1d into bby to finish:
  bby(0,ix)=zero
  do iy=1,nym1
    bby(iy,ix)=by1d(iy)
  enddo
  bby(ny,ix)=zero

enddo !End loop over columns

return
end subroutine

!==================================================

end module     
