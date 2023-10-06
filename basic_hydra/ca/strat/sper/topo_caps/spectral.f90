module spectral

use constants
use sta2dfft
use deriv1d

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Maximum x wavenumber:
integer,parameter:: nwx=nx/2,nwxm1=nwx-1,nwxp1=nwx+1

 !Common arrays, constants:
double precision:: pbar(0:ny)
double precision:: decy(nym1,0:nxm1)
double precision:: yh0(0:ny),yh1(0:ny)
double precision:: green(0:nxm1,0:ny),diss(0:nxm1,0:ny)
double precision:: rkx(0:nxm1),hrkx(nx),rky(ny)
double precision:: dafx(0:nx),dafy(0:ny) 
double precision:: frkx(0:nxm1),fhrkx(nx),frky(ny)
double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

 !Arrays related to the conformal map:
double precision:: yori(0:ny,0:nxm1),xori(0:ny,0:nxm1)
double precision:: dyoridx(0:ny,0:nxm1),dyoridy(0:ny,0:nxm1)
double precision:: confac(0:ny,0:nxm1),confaci(0:ny,0:nxm1)

!====================================================================!
! From main code: call init_invert            to initialise          !
! then            call main_invert(zz,uu,vv)  to perform inversion   !
!====================================================================!

contains

!===========================
subroutine init_spectral(bbdif,hh)

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: rkx1d(0:nxm1),dhdx(0:nxm1),hh(0:nxm1)

!----------------------

 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,rky)
call init_deriv(nx,ellx,rkx1d)

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
do kx=1,nxm1
  wratx=hrkx(kx)/rkxmax
  fhrkx(kx)=hrkx(kx)*exp(-36.d0*wratx**36.d0)
enddo

 !Define y wavenumbers:
scy=pi/elly
rkymax=scy*dble(ny)
do ky=1,ny
  wraty=rky(ky)/rkymax
  frky(ky)=rky(ky)*exp(-36.d0*wraty**36.d0)
enddo
 
 !Define de-aliasing filter (2/3 rule):
dafx(0)=one
do kx=1,nxm1
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
do kx=1,nxm1
  diss(kx,0)=f12*visc*rkx(kx)**(2*nnu)
enddo
do ky=1,ny
  diss(0,ky)=f12*visc*rky(ky)**(2*nnu)
enddo
do ky=1,ny
  do kx=1,nxm1
    diss(kx,ky)=f12*visc*(rkx(kx)**2+rky(ky)**2)**nnu
  enddo
enddo

 !Form arrays for conformal mapping:
call forfft(1,nx,hh,xtrig,xfactors)
call deriv(1,nx,rkx1d,hh,dhdx)
 !Define the interior semi-spectral field:
do iy=0,ny
  xori(iy,0)=zero
  yori(iy,0)=hh(0)*yh0(iy)
  dyoridx(iy,0)=zero
  dyoridy(iy,0)=zero
enddo
do kx=1,nxm1
  fac=rkx(kx)*elly
  divy=one/sinh(fac)
  divx=divy/rkx(kx)
  divd=-divy*rkx(kx)
  do iy=0,ny
    xori(iy,kx)=   dhdx(kx)*cosh(fac*yh0(iy))*divx
    yori(iy,kx)=     hh(kx)*sinh(fac*yh0(iy))*divy
    dyoridx(iy,kx)=dhdx(kx)*sinh(fac*yh0(iy))*divy
    dyoridy(iy,kx)=  hh(kx)*cosh(fac*yh0(iy))*divd
  enddo
enddo
 !Invert using a full transform in x:
call revfft(ny+1,nx,xori,xtrig,xfactors)
call revfft(ny+1,nx,yori,xtrig,xfactors)
call revfft(ny+1,nx,dyoridx,xtrig,xfactors)
call revfft(ny+1,nx,dyoridy,xtrig,xfactors)

 !Add on the x,y coordinate to the deviations:
do ix=0,nxm1
  do iy=0,ny
    xori(iy,ix)=xori(iy,ix)+glx*dble(ix)+xmin
    yori(iy,ix)=yori(iy,ix)+gly*dble(iy)+ymin
    dyoridy(iy,ix)=dyoridy(iy,ix)+one
  enddo
enddo

 !Define conformal factor for Poisson problem, and rescale derivatives 
 !of transform for use in calculating grad(b) in conformal space:
do ix=0,nxm1
  do iy=0,ny
    confac(iy,ix)=dyoridx(iy,ix)**2+dyoridy(iy,ix)**2
    confaci(iy,ix)=one/confac(iy,ix)
    dyoridx(iy,ix)=dyoridx(iy,ix)*confaci(iy,ix)
    dyoridy(iy,ix)=dyoridy(iy,ix)*confaci(iy,ix)
  enddo
enddo

return 
end subroutine

!===========================
subroutine main_invert(zz,uu,vv,t)

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: zz(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
 !Local arrays:
double precision:: ss(0:nxm1,0:ny),pp(0:ny,0:nxm1)
double precision:: pbot(0:nxm1),ptop(0:nxm1),cppy(nym1,0:nxm1)

!----------------------------------------------------------
 !Apply conformal factor to vorticity (use pp temporarily):
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix)=zz(iy,ix)*confac(iy,ix)
  enddo
enddo

 !Solve for psi (pp):

 !(1) compute mean vorticity (zbar):
zbar=zero
do ix=0,nxm1
  zbar=zbar+pp(0,ix)+pp(ny,ix)
enddo
zbar=f12*zbar
do ix=0,nxm1
  do iy=1,nym1
    zbar=zbar+pp(iy,ix)
  enddo
enddo
zbar=zbar/dble(nx*ny)

 !(2) Remove mean vorticity from pp:
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix)=pp(iy,ix)-zbar
  enddo
enddo
 
 !(3) FFT pp (vorticity) and invert to get uncorrected streamfunction pp:
call ptospc_fc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)
do ky=0,ny
  do kx=0,nxm1
    ss(kx,ky)=green(kx,ky)*ss(kx,ky)
  enddo
enddo
call spctop_fc(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

 !(4) Add part of pp (streamfunction) due to mean vorticity:
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

call getvel(pp,uu,vv,t)

return
end subroutine

!===================================================================

subroutine getvel(pp,uu,vv,t)
! Computes the velocity components uu & vv from the streamfunction
! pp via uu = -lambda^{-1}d(pp)/dy and vv = lambda^{-1}d(pp)/dx,
! where lambda is the conformal factor (confac).
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

 !Apply de-aliasing filter:
do ky=1,ny
  do kx=0,nxm1
    pps(kx,ky)=pps(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo

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

 !Copy -uui into uu, add mean reference velocity and tide, 
 !then apply conformal factor:
uadd=uref+utidemax*sin(ftide*t)
do ix=0,nxm1
  do iy=0,ny
    uu(iy,ix)=confaci(iy,ix)*(uadd-uui(iy,ix))
    vv(iy,ix)=confaci(iy,ix)*vv(iy,ix)
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
