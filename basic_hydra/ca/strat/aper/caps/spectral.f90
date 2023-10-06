module spectral

use constants
use sta2dfft

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Common arrays, constants:
double precision:: pbar(0:ny,0:nx)
double precision:: decx(nxm1,nym1),decy(nym1,nxm1)
double precision:: xh0(0:nx),xh1(0:nx)
double precision:: yh0(0:ny),yh1(0:ny)
double precision:: green(0:nx,0:ny),diss(0:nx,0:ny)
double precision:: dafx(0:nx),dafy(0:ny) 

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

double precision:: rkx(nx),rky(ny)
double precision:: rkxf(nx),rkyf(ny)

double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(0:nx,0:ny),kmax,kmaxred

!====================================================================!
! From main code: call init_invert            to initialise          !
! then            call main_invert(zz,uu,vv)  to perform inversion   !
!====================================================================!

contains

!===========================
subroutine init_spectral(bbdif)

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

!----------------------

 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,rkx,rky)

 !Weights for removing corner values of the streamfunction:
fac=one/dble(nx)
do ix=0,nx
  xh1(ix)=fac*dble(ix)
  xh0(ix)=one-xh1(ix)
enddo

fac=one/dble(ny)
do iy=0,ny
  yh1(iy)=fac*dble(iy)
  yh0(iy)=one-yh1(iy)
enddo

 !Define part of streamfunction proportional to the mean vorticity:
do ix=0,nx
  do iy=0,ny
    pbar(iy,ix)=-f14*(ellx**2*xh0(ix)*xh1(ix)+elly**2*yh0(iy)*yh1(iy))
  enddo
enddo

 !Define x wavenumbers:
scx=pi/ellx
rkxmax=scx*dble(nx)
do kx=1,nx
  rkx(kx)=scx*dble(kx)
  wratx=rkx(kx)/rkxmax
  rkxf(kx)=rkx(kx)*exp(-36.d0*wratx**36.d0)
enddo

 !Define y wavenumbers:
scy=pi/elly
rkymax=scy*dble(ny)
do ky=1,ny
  rky(ky)=scy*dble(ky)
  wraty=rky(ky)/rkymax
  rkyf(ky)=rky(ky)*exp(-36.d0*wraty**36.d0)
enddo

 !Define de-aliasing filter (2/3 rule):
dafx(0)=one
do kx=1,nx
  if (rkx(kx) .lt. f23*rkxmax) then
    dafx(kx)=one
  else
    dafx(kx)=zero
  endif
enddo

dafy(0)=one
do ky=1,ny
  if (rky(ky) .lt. f23*rkymax) then
    dafy(ky)=one
  else
    dafy(ky)=zero
  endif
enddo
 
 !Initialise arrays for computing the spectrum of any field:
delk=sqrt(scx**2+scy**2)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do ky=0,ny
  do kx=0,nx
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
do kx=1,nx
  green(kx,0)=-one/rkx(kx)**2
enddo
do ky=1,ny
  green(0,ky)=-one/rky(ky)**2
enddo
do ky=1,ny
  do kx=1,nx
    green(kx,ky)=-one/(rkx(kx)**2+rky(ky)**2)
  enddo
enddo

 !Hyperbolic functions used for solutions of Laplace's equation:
do kx=1,nxm1
  fac=rkx(kx)*elly
  div=one/(one-exp(-two*fac))
  do iy=1,nym1
    argm=fac*(one-yh1(iy))
    argp=fac*(one+yh1(iy))
    decy(iy,kx)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

do ky=1,nym1
  fac=rky(ky)*ellx
  div=one/(one-exp(-two*fac))
  do ix=1,nxm1
    argm=fac*(one-xh1(ix))
    argp=fac*(one+xh1(ix))
    decx(ix,ky)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

 !Define dissipation operator:
if (nnu .eq. 1) then 
  visc=prediss*sqrt(bbdif)/(rkxmax**2+rkymax**2)**0.75d0
   !bbdif=bb_max-bb_min
else
  visc=prediss/(rkxmax**2+rkymax**2)**nnu
endif

 !Write viscosity to log file:
write(*,*)
write(*,'(a,1x,1p,e14.7)') ' viscosity = ',visc
write(*,*)

diss(0,0)=zero
do kx=1,nx
  diss(kx,0)=f12*visc*rkx(kx)**(2*nnu)
enddo
do ky=1,ny
  diss(0,ky)=f12*visc*rky(ky)**(2*nnu)
enddo
do ky=1,ny
  do kx=1,nx
    diss(kx,ky)=f12*visc*(rkx(kx)**2+rky(ky)**2)**nnu
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
double precision:: zz(0:ny,0:nx),uu(0:ny,0:nx),vv(0:ny,0:nx)
 !Local arrays:
double precision:: ss(0:nx,0:ny),pp(0:ny,0:nx)
double precision:: pbot(nx),ptop(nx),cppy(nym1,nx)
double precision:: plft(ny),prgt(ny),cppx(nxm1,ny)

!--------------------------------
 !Solve for psi (pp):

 !(1) compute mean zz (zbar):
zbar=zero
do ix=1,nxm1
  zbar=zbar+zz(0,ix)+zz(ny,ix)
enddo
do iy=1,nym1
  zbar=zbar+zz(iy,0)+zz(iy,nx)
enddo
zbar=f12*zbar+f14*(zz(0,0)+zz(ny,0)+zz(0,nx)+zz(ny,nx))
do ix=1,nxm1
  do iy=1,nym1
    zbar=zbar+zz(iy,ix)
  enddo
enddo
zbar=zbar/dble(nx*ny)

 !(2) Remove mean vorticity from zz:
do ix=0,nx
  do iy=0,ny
    pp(iy,ix)=zz(iy,ix)-zbar
  enddo
enddo
 
 !(3) FFT zz and invert to get uncorrected streamfunction pp:
call ptospc_cc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)
do ky=0,ny
  do kx=0,nx
    ss(kx,ky)=green(kx,ky)*ss(kx,ky)
  enddo
enddo
call spctop_cc(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

 !(4) Add part of pp due to mean vorticity:
do ix=0,nx
  do iy=0,ny
    pp(iy,ix)=pp(iy,ix)+zbar*pbar(iy,ix)
  enddo
enddo

 !(5) Remove a bi-linear function so that pp is zero at the corners:
sw00=pp(0,0)
sw10=pp(ny,0)
sw01=pp(0,nx)
sw11=pp(ny,nx)
do ix=0,nx
  do iy=0,ny
    pp(iy,ix)=pp(iy,ix)-(sw00*xh0(ix)+sw01*xh1(ix))*yh0(iy)&
                      &-(sw10*xh0(ix)+sw11*xh1(ix))*yh1(iy)
  enddo
enddo
 !Note:  xh0 = (xmax - x)/ellx, xh1 = (x - xmin)/ellx etc.

 !(6) Do a sine transform of pp at y = ymin and ymax and obtain the
 !    interior field (cppy) that must be subtracted to give pp = 0
 !    at y = ymin and ymax:
do ix=1,nxm1
  pbot(ix)=pp(0,ix)
  ptop(ix)=pp(ny,ix)
enddo
call dst(1,nx,pbot,xtrig,xfactors)
call dst(1,nx,ptop,xtrig,xfactors)

 !Define the interior semi-spectral field:
do kx=1,nxm1
  do iy=1,nym1
    cppy(iy,kx)=pbot(kx)*decy(ny-iy,kx)+ptop(kx)*decy(iy,kx)
  enddo
enddo
 !Invert using a sine transform:
call dst(nym1,nx,cppy,xtrig,xfactors)

 !(7) Do a sine transform of pp at x = xmin and xmax and obtain the
 !    interior field (cppx) that must be subtracted to give pp = 0
 !    at x = xmin and xmax:
do iy=1,nym1
  plft(iy)=pp(iy,0)
  prgt(iy)=pp(iy,nx)
enddo
call dst(1,ny,plft,ytrig,yfactors)
call dst(1,ny,prgt,ytrig,yfactors)

 !Define the interior semi-spectral field:
do ky=1,nym1
  do ix=1,nxm1
    cppx(ix,ky)=plft(ky)*decx(nx-ix,ky)+prgt(ky)*decx(ix,ky)
  enddo
enddo
 !Invert using a sine transform:
call dst(nxm1,ny,cppx,ytrig,yfactors)

 !(8) Remove cppx and cppy to obtain the final streamfunction pp:
do iy=0,ny
  pp(iy,0 )=zero
  pp(iy,nx)=zero
enddo
do ix=0,nx
  pp(0, ix)=zero
  pp(ny,ix)=zero
enddo

do ix=1,nxm1
  do iy=1,nym1
    pp(iy,ix)=pp(iy,ix)-cppx(ix,iy)-cppy(iy,ix)
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
double precision:: pp(0:ny,0:nx),uu(0:ny,0:nx),vv(0:ny,0:nx)

 !Local arrays:
double precision:: ppi(ny,nx),pps(nx,ny)
double precision:: ppx(0:nx,ny),vvi(ny,0:nx)
double precision:: ppy(nx,0:ny),uui(0:ny,nx)

 !Copy non-zero interior values of pp to ppi:
do ix=1,nxm1
  do iy=1,nym1
    ppi(iy,ix)=pp(iy,ix)
  enddo
enddo

 !Transform ppi to spectral space:
call ptospc_ss(nx,ny,ppi,pps,xfactors,yfactors,xtrig,ytrig)

 !Apply de-aliasing filter:
do ky=1,ny
  do kx=1,nx
    pps(kx,ky)=pps(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo

 !Compute d(ppi)/dx = ppx spectrally:
call xderiv_ss(nx,ny,rkx,pps,ppx)

 !Transform ppx back to physical space as vvi:
call spctop_cs(nx,ny,ppx,vvi,xfactors,yfactors,xtrig,ytrig)

 !Copy vvi into vv and add on zero edge values at iy = 0 & ny:
do ix=0,nx
  vv(0,ix)=zero
  do iy=1,nym1
    vv(iy,ix)=vvi(iy,ix)
  enddo
  vv(ny,ix)=zero
enddo

 !Compute d(ppi)/dy = ppy spectrally:
call yderiv_ss(nx,ny,rky,pps,ppy)

 !Transform ppy back to physical space as uui:
call spctop_sc(nx,ny,ppy,uui,xfactors,yfactors,xtrig,ytrig)

 !Copy -uui into uu and add on zero edge values at ix = 0 & nx:
do ix=1,nxm1
  do iy=0,ny
    uu(iy,ix)=-uui(iy,ix)
  enddo
enddo
do iy=0,ny
  uu(iy, 0)=zero
  uu(iy,nx)=zero
enddo

return
end subroutine

!===================================================================

subroutine spec1d_cc(rvar,spec)
! Computes the 1d spectrum of a field rvar which is
! represented by a cosine series in x and a cosine series in y
! and returns the result in spec.
!   *** Warning: rvar is modified by this routine.***

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(0:ny,0:nx),spec(0:max(nx,ny))
 !Local array:
double precision:: wks(0:nx,0:ny)

 !Transform rvar to spectral space:
call ptospc_cc(nx,ny,rvar,wks,xfactors,yfactors,xtrig,ytrig)

do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(0,0)
spec(k)=spec(k)+f14*wks(0,0)**2

 !y-independent mode:
do kx=1,nx
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
  do kx=1,nx
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
double precision:: bb(0:ny,0:nx),bbx(0:ny,0:nx)

 !Flag needed for FFTs:
ifail=0

 !Copy bb into bbx for FFT:
do ix=0,nx
  do iy=0,ny
    bbx(iy,ix)=bb(iy,ix)
  enddo
enddo

 !Carry out an x cosine transform on bbx:
call dct(nyp1,nx,bbx,xtrig,xfactors)

 !Take derivative spectrally by wavenumber multiplication:
do kx=1,nxm1
  do iy=0,ny
    bbx(iy,kx)=-rkxf(kx)*bbx(iy,kx)
  enddo
enddo
do iy=0,ny
  bbx(iy,nx)=zero
enddo

 !Carry out an x sine transform on bbx:
call dst(nyp1,nx,bbx(0,1),xtrig,xfactors)

 !Add zero edge values:
do iy=0,ny
  bbx(iy, 0)=zero
  bbx(iy,nx)=zero
enddo

return
end subroutine

!============================================================
subroutine byderiv(bb,bby)

! Subroutine to take the y derivative of a field bb, assumed to
! be even across the y boundaries.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: bb(0:ny,0:nx),bby(0:ny,0:nx)

 !Local array:
double precision:: b1d(0:ny),by1d(ny)

 !Flag needed for FFTs:
ifail=0

 !Loop over columns, computing db/dy along each x grid line:
do ix=0,nx
   !Copy bb into b1d for FFT:
  do iy=0,ny
    b1d(iy)=bb(iy,ix)
  enddo

   !Carry out a y cosine transform on b1d:
  call dct(1,ny,b1d,ytrig,yfactors)

   !Take derivative spectrally by wavenumber multiplication:
  do ky=1,nym1
    by1d(ky)=-rkyf(ky)*b1d(ky)
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

enddo !End loop over rows

return
end subroutine

!============================================================
end module     
