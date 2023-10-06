module spectral

use constants

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Common arrays, constants:
double precision:: pbar(0:ny,0:nx)
double precision:: decx(nxm1,nym1),decy(nym1,nxm1)
double precision:: xh0(0:nx),xh1(0:nx)
double precision:: yh0(0:ny),yh1(0:ny)
double precision:: green(0:nx,0:ny)
double precision:: cosxtab(2*nx),cosytab(2*ny)
double precision:: sinxtab(2*nx),sinytab(2*ny)
double precision:: wk1d(0:nx+ny-2),wk2d(0:ny,0:nx)
double precision:: rkx(nx),rky(ny)
double precision:: rkxf(nx),rkyf(ny)
double precision:: spmf(0:(nx+ny)/2),alk((nx+ny)/2)
integer:: kmag(0:nx,0:ny),kmax,kmaxred,ifail

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

 !Set up FFT tables:
ifail=0
call c06hbf(nyp1,nx,pbar, 'i',cosxtab,wk2d,ifail)
call c06hbf(nxp1,ny,green,'i',cosytab,wk2d,ifail)
call c06haf(   1,nx,rkx,  'i',sinxtab,wk1d,ifail)
call c06haf(   1,ny,rky,  'i',sinytab,wk1d,ifail)
 !Note - pbar,..,rky are used as dummy arrays 
 !     - they are overwritten below 

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
 
 !Initialise arrays for computing the spectrum of any field:
delk=sqrt(scx*scy)
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
call ptospc_cc(pp,ss)
do ky=0,ny
  do kx=0,nx
    ss(kx,ky)=green(kx,ky)*ss(kx,ky)
  enddo
enddo
call spctop_cc(ss,pp)

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
call c06haf(1,nx,pbot,'s',sinxtab,wk1d,ifail)
call c06haf(1,nx,ptop,'s',sinxtab,wk1d,ifail)

 !Define the interior semi-spectral field:
do kx=1,nxm1
  do iy=1,nym1
    cppy(iy,kx)=pbot(kx)*decy(ny-iy,kx)+ptop(kx)*decy(iy,kx)
  enddo
enddo
 !Invert using a sine transform:
call c06haf(nym1,nx,cppy,'s',sinxtab,wk2d,ifail)

 !(7) Do a sine transform of pp at x = xmin and xmax and obtain the
 !    interior field (cppx) that must be subtracted to give pp = 0
 !    at x = xmin and xmax:
do iy=1,nym1
  plft(iy)=pp(iy,0)
  prgt(iy)=pp(iy,nx)
enddo
call c06haf(1,ny,plft,'s',sinytab,wk1d,ifail)
call c06haf(1,ny,prgt,'s',sinytab,wk1d,ifail)

 !Define the interior semi-spectral field:
do ky=1,nym1
  do ix=1,nxm1
    cppx(ix,ky)=plft(ky)*decx(nx-ix,ky)+prgt(ky)*decx(ix,ky)
  enddo
enddo
 !Invert using a sine transform:
call c06haf(nxm1,ny,cppx,'s',sinytab,wk2d,ifail)

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
call ptospc_ss(ppi,pps)

 !Compute d(ppi)/dx = ppx spectrally:
call xderiv_ss(pps,ppx)

 !Transform ppx back to physical space as vvi:
call spctop_cs(ppx,vvi)

 !Copy vvi into vv and add on zero edge values at iy = 0 & ny:
do ix=0,nx
  vv(0,ix)=zero
  do iy=1,nym1
    vv(iy,ix)=vvi(iy,ix)
  enddo
  vv(ny,ix)=zero
enddo

 !Compute d(ppi)/dy = ppy spectrally:
call yderiv_ss(pps,ppy)

 !Transform ppy back to physical space as uui:
call spctop_sc(ppy,uui)

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
double precision:: rvar(0:ny,0:nx),spec(0:(nx+ny)/2)
 !Local array:
double precision:: wks(0:nx,0:ny)

 !Transform rvar to spectral space:
call ptospc_cc(rvar,wks)

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

 !Copy bb into bbx for FFT:
do ix=0,nx
  do iy=0,ny
    bbx(iy,ix)=bb(iy,ix)
  enddo
enddo

 !Carry out an x cosine transform on bbx:
call c06hbf(nyp1,nx,bbx,'s',cosxtab,wk2d,ifail)

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
call c06haf(nyp1,nx,bbx(0,1),'s',sinxtab,wk2d,ifail)  

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

 !Loop over columns, computing db/dy along each x grid line:
do ix=0,nx
   !Copy bb into b1d for FFT:
  do iy=0,ny
    b1d(iy)=bb(iy,ix)
  enddo

   !Carry out a y cosine transform on b1d:
  call c06hbf(1,ny,b1d,'s',cosytab,wk1d,ifail)

   !Take derivative spectrally by wavenumber multiplication:
  do ky=1,nym1
    by1d(ky)=-rkyf(ky)*b1d(ky)
  enddo
  by1d(ny)=zero

   !Carry out a y sine transform on by1d:
  call c06haf(1,ny,by1d,'s',sinytab,wk1d,ifail)

   !Copy by1d into bby to finish:
  bby(0,ix)=zero
  do iy=1,nym1
    bby(iy,ix)=by1d(iy)
  enddo
  bby(ny,ix)=zero

enddo !End loop over rows

return
end subroutine

!==================================================

subroutine ptospc_cc(rvar,svar)
! Performs a physical -> spectral transform of a variable rvar
! represented by a cosine series in x and a cosine series in y
! and returns the result (transposed) in svar (with kx as the
! inner dimension).

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(0:ny,0:nx),svar(0:nx,0:ny)

 !Carry out x cosine transform first:
call c06hbf(nyp1,nx,rvar,'s',cosxtab,wk2d,ifail)

 !Transpose array:
do kx=0,nx
  do iy=0,ny
    svar(kx,iy)=rvar(iy,kx)
  enddo
enddo

 !Carry out y cosine transform on transposed array:
call c06hbf(nxp1,ny,svar,'s',cosytab,wk2d,ifail)

return
end subroutine

!===================================================================

subroutine ptospc_cs(rvar,svar)
! Performs a physical -> spectral transform of a variable rvar
! represented by a cosine series in x and a  sine  series in y
! and returns the result (transposed) in svar (with kx as the
! inner dimension).

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(ny,0:nx),svar(0:nx,ny)

 !Carry out x cosine transform first:
call c06hbf(ny,nx,rvar,'s',cosxtab,wk2d,ifail)

 !Transpose array:
do kx=0,nx
  do iy=1,nym1
    svar(kx,iy)=rvar(iy,kx)
  enddo
  svar(kx,ny)=zero
enddo

 !Carry out y sine transform on transposed array:
call c06haf(nxp1,ny,svar,'s',sinytab,wk2d,ifail)

return
end subroutine

!===================================================================

subroutine ptospc_sc(rvar,svar)
! Performs a physical -> spectral transform of a variable rvar
! represented by a  sine  series in x and a cosine series in y
! and returns the result (transposed) in svar (with kx as the
! inner dimension).

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(0:ny,nx),svar(nx,0:ny)

 !Carry out x sine transform first:
call c06haf(nyp1,nx,rvar,'s',sinxtab,wk2d,ifail)

 !Transpose array:
do kx=1,nxm1
  do iy=0,ny
    svar(kx,iy)=rvar(iy,kx)
  enddo
enddo
do iy=0,ny
  svar(nx,iy)=zero
enddo

 !Carry out y cosine transform on transposed array:
call c06hbf(nx,ny,svar,'s',cosytab,wk2d,ifail)

return
end subroutine

!===================================================================

subroutine ptospc_ss(rvar,svar)
! Performs a physical -> spectral transform of a variable rvar
! represented by a  sine  series in x and a  sine  series in y
! and returns the result (transposed) in svar (with kx as the
! inner dimension).

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(ny,nx),svar(nx,ny)

 !Carry out x sine transform first:
call c06haf(ny,nx,rvar,'s',sinxtab,wk2d,ifail)

 !Transpose array:
do kx=1,nxm1
  do iy=1,nym1
    svar(kx,iy)=rvar(iy,kx)
  enddo
  svar(kx,ny)=zero
enddo
do iy=1,ny
  svar(nx,iy)=zero
enddo

 !Carry out y sine transform on transposed array:
call c06haf(nx,ny,svar,'s',sinytab,wk2d,ifail)

return
end subroutine

!===================================================================

subroutine spctop_cc(svar,rvar)
! Performs a spectral -> physical transform of a variable svar
! represented by a cosine series in x and a cosine series in y
! and returns the result (transposed) in rvar (with iy as the
! inner dimension).

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(0:ny,0:nx),svar(0:nx,0:ny)

 !Carry out y cosine transform first:
call c06hbf(nxp1,ny,svar,'s',cosytab,wk2d,ifail)

 !Transpose array:
do kx=0,nx
  do iy=0,ny
    rvar(iy,kx)=svar(kx,iy)
  enddo
enddo

 !Carry out x cosine transform on transposed array:
call c06hbf(nyp1,nx,rvar,'s',cosxtab,wk2d,ifail)

return
end subroutine

!===================================================================

subroutine spctop_cs(svar,rvar)
! Performs a spectral -> physical transform of a variable svar
! represented by a cosine series in x and a  sine  series in y
! and returns the result (transposed) in rvar (with iy as the
! inner dimension).

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(ny,0:nx),svar(0:nx,ny)

 !Carry out y sine transform first:
call c06haf(nxp1,ny,svar,'s',sinytab,wk2d,ifail)

 !Transpose array:
do kx=0,nx
  do iy=1,nym1
    rvar(iy,kx)=svar(kx,iy)
  enddo
  rvar(ny,kx)=zero
enddo

 !Carry out x cosine transform on transposed array:
call c06hbf(ny,nx,rvar,'s',cosxtab,wk2d,ifail)

return
end subroutine

!===================================================================

subroutine spctop_sc(svar,rvar)
! Performs a spectral -> physical transform of a variable svar
! represented by a  sine  series in x and a cosine series in y
! and returns the result (transposed) in rvar (with iy as the
! inner dimension).

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(0:ny,nx),svar(nx,0:ny)

 !Carry out y cosine transform first:
call c06hbf(nx,ny,svar,'s',cosytab,wk2d,ifail)

 !Transpose array:
do kx=1,nxm1
  do iy=0,ny
    rvar(iy,kx)=svar(kx,iy)
  enddo
enddo
do iy=0,ny
  rvar(iy,nx)=zero
enddo

 !Carry out x sine transform on transposed array:
call c06haf(nyp1,nx,rvar,'s',sinxtab,wk2d,ifail)

return
end subroutine

!===================================================================

subroutine spctop_ss(svar,rvar)
! Performs a spectral -> physical transform of a variable svar
! represented by a  sine  series in x and a  sine  series in y
! and returns the result (transposed) in rvar (with iy as the
! inner dimension).

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(ny,nx),svar(nx,ny)

 !Carry out y sine transform first:
call c06haf(nx,ny,svar,'s',sinytab,wk2d,ifail)

 !Transpose array:
do kx=1,nxm1
  do iy=1,nym1
    rvar(iy,kx)=svar(kx,iy)
  enddo
  rvar(ny,kx)=zero
enddo
do iy=1,ny
  rvar(iy,nx)=zero
enddo

 !Carry out x sine transform on transposed array:
call c06haf(ny,nx,rvar,'s',sinxtab,wk2d,ifail)

return
end subroutine

!===================================================================

subroutine xderiv_ss(var,der)
! Computes der = d(var)/dx, spectrally, for a variable var
! represented by a  sine  series in x and a  sine  series in y.
! *** both var and der are spectral ***
! ==> der changes to a cosine series in x
! @@@ the mean value of der (0 wavenumber) is assumed to be 0

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(nx,ny),der(0:nx,ny)

 !Carry out differentiation by wavenumber multiplication:
do ky=1,nym1
  der(0 ,ky)=zero
   !The above implies the mean value of der for every ky is zero.
  do kx=1,nxm1
    der(kx,ky)=rkx(kx)*var(kx,ky)
  enddo
  der(nx,ky)=zero
   !der = 0 when kx = nx since var = 0 when kx = nx.
enddo

return
end subroutine

!===================================================================

subroutine yderiv_ss(var,der)
! Computes der = d(var)/dy, spectrally, for a variable var
! represented by a  sine  series in x and a  sine  series in y.
! *** both var and der are spectral ***
! ==> der changes to a cosine series in y
! @@@ the mean value of der (0 wavenumber) is assumed to be 0

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(nx,ny),der(nx,0:ny)

 !Carry out differentiation by wavenumber multiplication:
do kx=1,nxm1
  der(kx,0)=zero
   !The above implies the mean value of der for every kx is zero.
  der(kx,ny)=zero
   !der = 0 when ky = ny since var = 0 when ky = ny.
enddo
do ky=1,nym1
  do kx=1,nxm1
    der(kx,ky)=rky(ky)*var(kx,ky)
  enddo
enddo

return
end subroutine

!===================================================================

subroutine xderiv_cs(var,der)
! Computes der = d(var)/dx, spectrally, for a variable var
! represented by a cosine series in x and a  sine  series in y.
! *** both var and der are spectral ***
! ==> der changes to a sine series in x

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(0:nx,ny),der(nx,ny)

 !Carry out differentiation by wavenumber multiplication:
do ky=1,nym1
  do kx=1,nxm1
    der(kx,ky)=-rkx(kx)*var(kx,ky)
  enddo
  der(nx,ky)=zero
enddo

return
end subroutine

!===================================================================

subroutine yderiv_cs(var,der)
! Computes der = d(var)/dy, spectrally, for a variable var
! represented by a cosine series in x and a  sine  series in y.
! *** both var and der are spectral ***
! ==> der changes to a cosine series in y
! @@@ the mean value of der (0 wavenumber) is assumed to be 0

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(0:nx,ny),der(0:nx,0:ny)

 !Carry out differentiation by wavenumber multiplication:
do kx=0,nx
  der(kx,0)=zero
   !The above implies the mean value of der for every kx is zero.
  der(kx,ny)=zero
   !der = 0 when ky = ny since var = 0 when ky = ny.
enddo
do ky=1,nym1
  do kx=0,nx
    der(kx,ky)=rky(ky)*var(kx,ky)
  enddo
enddo

return
end subroutine

!===================================================================

subroutine xderiv_sc(var,der)
! Computes der = d(var)/dx, spectrally, for a variable var
! represented by a  sine  series in x and a cosine series in y.
! *** both var and der are spectral ***
! ==> der changes to a cosine series in x
! @@@ the mean value of der (0 wavenumber) is assumed to be 0

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(nx,0:ny),der(0:nx,0:ny)

 !Carry out differentiation by wavenumber multiplication:
do ky=0,ny
  der(0 ,ky)=zero
   !The above implies the mean value of der for every ky is zero.
  do kx=1,nxm1
    der(kx,ky)=rkx(kx)*var(kx,ky)
  enddo
  der(nx,ky)=zero
   !der = 0 when kx = nx since var = 0 when kx = nx.
enddo

return
end subroutine

!===================================================================

subroutine yderiv_sc(var,der)
! Computes der = d(var)/dy, spectrally, for a variable var
! represented by a  sine  series in x and a cosine series in y.
! *** both var and der are spectral ***
! ==> der changes to a sine series in y

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(nx,0:ny),der(nx,ny)

 !Carry out differentiation by wavenumber multiplication:
do kx=1,nxm1
  der(kx,ny)=zero
enddo
do ky=1,nym1
  do kx=1,nxm1
    der(kx,ky)=-rky(ky)*var(kx,ky)
  enddo
enddo

return
end subroutine

!===================================================================

subroutine xderiv_cc(var,der)
! Computes der = d(var)/dx, spectrally, for a variable var
! represented by a cosine series in x and a cosine series in y.
! *** both var and der are spectral ***
! ==> der changes to a sine series in x

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(0:nx,0:ny),der(nx,0:ny)

 !Carry out differentiation by wavenumber multiplication:
do ky=0,ny
  do kx=1,nxm1
    der(kx,ky)=-rkx(kx)*var(kx,ky)
  enddo
  der(nx,ky)=zero
enddo

return
end subroutine

!===================================================================

subroutine yderiv_cc(var,der)
! Computes der = d(var)/dy, spectrally, for a variable var
! represented by a cosine series in x and a cosine series in y.
! *** both var and der are spectral ***
! ==> der changes to a sine series in y

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(0:nx,0:ny),der(0:nx,ny)

 !Carry out differentiation by wavenumber multiplication:
do kx=0,nx
  der(kx,ny)=zero
enddo
do ky=1,nym1
  do kx=0,nx
    der(kx,ky)=-rky(ky)*var(kx,ky)
  enddo
enddo

return
end subroutine

!===================================================================

end module     
