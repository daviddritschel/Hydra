module spectral

use constants
use generic
use sta2dfft

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Maximum x wavenumber:
integer,parameter:: nwx=nx/2,nwxm1=nwx-1,nwxp1=nwx+1

 !Common arrays, constants:
double precision:: yg(0:ny),yh0(0:ny),yh1(0:ny),bety(0:ny)
double precision:: decy1(nym1,nxm1),decy2(nym1,nxm1)
double precision:: qgop1(0:nxm1,0:ny),qgop2(0:nxm1,0:ny),laplace(0:nxm1,0:ny)
double precision:: rkx(0:nxm1),hrkx(nx),rky(ny)

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

double precision:: frkx(0:nxm1),frky(ny)
double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(0:nxm1,0:ny),kmax,kmaxred

!==========================================================================!
! From main code: call init_invert                     to initialise       !
! then            call main_invert(qq,fhb,t,uu,vv,pp)  to perform inversion!
!==========================================================================!

contains

!===========================
subroutine init_spectral

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: ld

!---------------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,rky)

 !Fractional y grid values: 
fac=one/dble(ny)
do iy=0,ny
  yh1(iy)=fac*dble(iy)
  yh0(iy)=one-yh1(iy)
enddo

 !Define y & beta*y:
do iy=0,ny
  yg(iy)=ymin+gly*dble(iy)
  bety(iy)=beta*yg(iy)
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
 
!----------------------------------------------------------------------
 !Initialise arrays for computing the spectrum of any field:
delk=sqrt(scx**2+scy**2)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do kx=0,nxm1
  k=nint(rkx(kx)*delki)
  kmag(kx,0)=k
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
snorm=four*pi/dble(nx*ny)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(delk*dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=nint(sqrt((rkxmax**2+rkymax**2)/two)*delki)

!----------------------------------------------------------------------
 !Define inverse spectral inversion operators:
if (barot) then
   !The first mode is barotropic (kd1 = 0):
  qgop1(0,0)=zero
else
  qgop1(0,0)=-one/(kd1sq+small)
endif
qgop2(0,0)=-one/kd2sq
laplace(0,0)=zero
do kx=1,nxm1
  rksq=rkx(kx)**2
  qgop1(kx,0)=-one/(rksq+kd1sq)
  qgop2(kx,0)=-one/(rksq+kd2sq)
  laplace(kx,0)=-rksq
enddo
do ky=1,ny
  rksq=rky(ky)**2
  qgop1(0,ky)=-one/(rksq+kd1sq)
  qgop2(0,ky)=-one/(rksq+kd2sq)
  laplace(0,ky)=-rksq
enddo
do ky=1,ny
  do kx=1,nxm1
    rksq=rkx(kx)**2+rky(ky)**2
    qgop1(kx,ky)=-one/(rksq+kd1sq)
    qgop2(kx,ky)=-one/(rksq+kd2sq)
    laplace(kx,ky)=-rksq
  enddo
enddo

!----------------------------------------------------------------------
 !Hyperbolic functions used to correct boundary conditions in inversion:
 !First mode:
do kx=1,nxm1
  fac=sqrt(rkx(kx)**2+kd1sq)*elly
  div=one/(one-exp(-two*fac))
  do iy=1,nym1
    argm=fac*(one-yh1(iy))
    argp=fac*(one+yh1(iy))
    decy1(iy,kx)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

 !Second mode (kd2 > kd1 is assumed):
do kx=1,nxm1
  fac=sqrt(rkx(kx)**2+kd2sq)*elly
  div=one/(one-exp(-two*fac))
  do iy=1,nym1
    argm=fac*(one-yh1(iy))
    argp=fac*(one+yh1(iy))
    decy2(iy,kx)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

return
end subroutine

!========================================
subroutine main_invert(qq,fhb,t,uu,vv,pp)

 !Input:  PV field (qq), scaled bottom topography f_0*H_b/H_1 (fhb) & time (t)
 !Output: Velocity field (uu,vv) and streamfunction field (pp)

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(0:ny,0:nxm1,2),pp(0:ny,0:nxm1,2)
double precision:: uu(0:ny,0:nxm1,2),vv(0:ny,0:nxm1,2)
double precision:: fhb(0:ny,0:nxm1)
 !Local arrays:
double precision:: cppy(nym1,0:nxm1)
double precision:: p1(0:ny,0:nxm1),p2(0:ny,0:nxm1)
double precision:: u1(0:ny,0:nxm1),u2(0:ny,0:nxm1)
double precision:: v1(  ny,0:nxm1),v2(  ny,0:nxm1)
double precision:: ss(0:nxm1,0:ny),pps(0:nxm1,ny),ppx(0:nxm1,ny)
double precision:: pbot(0:nxm1),ptop(0:nxm1)
double precision:: p1bar(0:ny),u1bar(0:ny),p2bar(0:ny),u2bar(0:ny)
double precision:: z1bar(0:ny),dpdy(ny)

!---------------------------------------------------------------------
 !Divide PV anomaly (q - beta*y) into modes:
do ix=0,nxm1
  do iy=0,ny
     !Add topography (f_0*H_b/H_1) to PV in lower layer:
    qq1=qq(iy,ix,1)-bety(iy)-fhb(iy,ix)
    qq2=qq(iy,ix,2)-bety(iy)
    p1(iy,ix)=vec11*qq1+vec12*qq2
    p2(iy,ix)=vec21*qq1+vec22*qq2
  enddo
enddo
 !Here, p1 and p2 are used temporarily for the modal PV anomaly fields

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !Solve for the mode 1 flow field (p1, u1, v1):

if (barot) then
   !Here, the first mode is barotropic (p1 contains the relative vorticity)

   !(1) compute zonal-mean vorticity (z1bar):
  do iy=0,ny
    z1bar(iy)=p1(iy,0)
  enddo
  do ix=1,nxm1
    do iy=0,ny
      z1bar(iy)=z1bar(iy)+p1(iy,ix)
    enddo
  enddo
  do iy=0,ny
    z1bar(iy)=z1bar(iy)*dnxi
  enddo
   !dnxi = 1/nx

   !(2) Integrate -z1bar to obtain the zonal-mean zonal velocity, u1bar:
  u1bar(0)=zero
  do iy=1,ny
    u1bar(iy)=u1bar(iy-1)-hgly*(z1bar(iy-1)+z1bar(iy))
  enddo
   !Remove mean (to enforce zero global momentum):
  u1avg=f12*u1bar(ny)
  do iy=1,nym1
    u1avg=u1avg+u1bar(iy)
  enddo
  u1avg=u1avg*dnyi
   !dnyi = 1/ny
  do iy=0,ny
    u1bar(iy)=u1bar(iy)-u1avg
  enddo

   !(3) Integrate -u1bar to obtain the zonal-mean streamfunction, p1bar:
  p1bar(0)=zero
  do iy=1,ny
    p1bar(iy)=p1bar(iy-1)-hgly*(u1bar(iy-1)+u1bar(iy))
  enddo
   !Remove mean (to enforce mass conservation):
  p1avg=f12*p1bar(ny)
  do iy=1,nym1
    p1avg=p1avg+p1bar(iy)
  enddo
  p1avg=p1avg*dnyi
   !dnyi = 1/ny
  do iy=0,ny
    p1bar(iy)=p1bar(iy)-p1avg
  enddo

   !(4) Remove zonal-mean vorticity from total vorticity, p1:
  do ix=0,nxm1
    do iy=0,ny
      p1(iy,ix)=p1(iy,ix)-z1bar(iy)
    enddo
  enddo
 
   !FFT p1 and invert to get uncorrected **non-zonal** streamfunction p1:
  call ptospc_fc(nx,ny,p1,ss,xfactors,yfactors,xtrig,ytrig)
  do ky=0,ny
     !The kx = 0 component has been removed above (so should be zero):
    ss(0,ky)=zero
    do kx=1,nxm1
      ss(kx,ky)=qgop1(kx,ky)*ss(kx,ky)
    enddo
  enddo
  call spctop_fc(nx,ny,ss,p1,xfactors,yfactors,xtrig,ytrig)

   !(5) Do a full transform of p1 at y = ymin and ymax and obtain the
   !    interior field (cppy) that must be subtracted to give p1 = 0
   !    at y = ymin and ymax:
  do ix=0,nxm1
    pbot(ix)=p1(0,ix)
    ptop(ix)=p1(ny,ix)
  enddo
  call forfft(1,nx,pbot,xtrig,xfactors)
  call forfft(1,nx,ptop,xtrig,xfactors)

   !Define the interior semi-spectral field:
  do iy=1,nym1
    cppy(iy,0)=zero
  enddo
  do kx=1,nxm1
    do iy=1,nym1
      cppy(iy,kx)=pbot(kx)*decy1(ny-iy,kx)+ptop(kx)*decy1(iy,kx)
    enddo
  enddo
   !Invert using a full transform in x:
  call revfft(nym1,nx,cppy,xtrig,xfactors)

   !(6) Remove cppy to obtain the final non-zonal streamfunction p1:
  do ix=0,nxm1
    p1(0, ix)=zero
    p1(ny,ix)=zero
  enddo

   !Here we also prepare for the next step:
  do ix=0,nxm1
    do iy=1,nym1
      p1(iy,ix)=p1(iy,ix)-cppy(iy,ix)
      v1(iy,ix)=p1(iy,ix)
    enddo
  enddo

   !(7) Compute barotropic **non-zonal** velocity field:

   !Transform p1 (in v1) to spectral space:
  call ptospc_fs(nx,ny,v1,pps,xfactors,yfactors,xtrig,ytrig)

   !Compute d(p1)/dx = ppx spectrally:
  call xderiv_fs(nx,ny,hrkx,pps,ppx)

   !Transform ppx back to physical space as v1:
  call spctop_fs(nx,ny,ppx,v1,xfactors,yfactors,xtrig,ytrig)

   !Compute d(p1)/dy = ss spectrally:
  call yderiv_fs(nx,ny,rky,pps,ss)

   !Transform ss back to physical space as u1:
  call spctop_fc(nx,ny,ss,u1,xfactors,yfactors,xtrig,ytrig)
   !Note, u1 here contains d(p1)/dy in physical space

   !(8) Correct the streamfunction and zonal velocity (add on zonal-mean):
  do ix=0,nxm1
    do iy=0,ny
      p1(iy,ix)=p1bar(iy)+p1(iy,ix)
      u1(iy,ix)=u1bar(iy)-u1(iy,ix)
    enddo
  enddo

else
   !Here, kd1 > 0; invert the mode 1 Helmholtz operator:

   !(1) remove kd1^2*u1hmean*y from mode 1 PV anomaly to enforce mean zonal
   !    velocity boundary conditions, u1 = u1hmeam at y = ymin & ymax:
  dq1dy=kd1sq*u1hmean
  do ix=0,nxm1
    do iy=0,ny
      p1(iy,ix)=p1(iy,ix)-dq1dy*yg(iy)
    enddo
  enddo

   !(2) FFT p1 and invert to get uncorrected streamfunction p1:
  call ptospc_fc(nx,ny,p1,ss,xfactors,yfactors,xtrig,ytrig)
  do ky=0,ny
    do kx=0,nxm1
      ss(kx,ky)=qgop1(kx,ky)*ss(kx,ky)
    enddo
  enddo

   !(3) Obtain zonal mean streamfunction and zonal velocity:
   !    Extract zonal mean part of p1 from ss and remove it from ss:
  sqdn=one/sqrt(dble(nx))
  do ky=0,ny
    p1bar(ky)=sqdn*ss(0,ky)
    ss(0,ky)=zero
  enddo
   !The sqrt(nx) factor is needed due to a change in the x transform.

   !Calculate a y derivative of the zonal part to get the zonal velocity,
   !apart from the uniform flow u1hmean (note dpdy = 0 on boundaries):
  call yderiv_fc(1,ny,rky,p1bar,dpdy)

   !Transform p1bar & dpdy back to physical space (cosine & sine transforms):
  call dct(1,ny,p1bar,ytrig,yfactors)
  call dst(1,ny, dpdy,ytrig,yfactors)

   !Correct zonal mean streamfunction:
  do iy=0,ny
    p1bar(iy)=p1bar(iy)-u1hmean*yg(iy)
  enddo

   !Correct zonal mean zonal velocity:
  u1bar(0)=u1hmean
  do iy=1,nym1
    u1bar(iy)=u1hmean-dpdy(iy)
  enddo
  u1bar(ny)=u1hmean

   !Transform non-zonal streamfunction back to physical space:
  call spctop_fc(nx,ny,ss,p1,xfactors,yfactors,xtrig,ytrig)

   !(4) Do a full transform of p1 at y = ymin and ymax and obtain the
   !    interior field (cppy) that must be subtracted to give p1 = 0
   !    at y = ymin and ymax:
  do ix=0,nxm1
    pbot(ix)=p1(0,ix)
    ptop(ix)=p1(ny,ix)
  enddo
  call forfft(1,nx,pbot,xtrig,xfactors)
  call forfft(1,nx,ptop,xtrig,xfactors)

   !Define the (non-zonal, kx > 0) interior semi-spectral field:
  do iy=1,nym1
    cppy(iy,0)=zero
  enddo
  do kx=1,nxm1
    do iy=1,nym1
      cppy(iy,kx)=pbot(kx)*decy1(ny-iy,kx)+ptop(kx)*decy1(iy,kx)
    enddo
  enddo
   !Invert using a full transform in x:
  call revfft(nym1,nx,cppy,xtrig,xfactors)

   !(5) Remove cppy to obtain the final mode 1 (non-zonal) streamfunction p1:
  do ix=0,nxm1
    p1(0, ix)=zero
    p1(ny,ix)=zero
  enddo

   !Here we also prepare for the next step:
  do ix=0,nxm1
    do iy=1,nym1
      p1(iy,ix)=p1(iy,ix)-cppy(iy,ix)
      v1(iy,ix)=p1(iy,ix)
    enddo
  enddo

   !(6) Compute mode 1 velocity field (correct for boundary conditions below):

   !Transform p1 (in v1) to spectral space:
  call ptospc_fs(nx,ny,v1,pps,xfactors,yfactors,xtrig,ytrig)

   !Compute d(p1)/dx = ppx spectrally:
  call xderiv_fs(nx,ny,hrkx,pps,ppx)

   !Transform ppx back to physical space as v1:
  call spctop_fs(nx,ny,ppx,v1,xfactors,yfactors,xtrig,ytrig)

   !Compute d(p1)/dy = ss spectrally:
  call yderiv_fs(nx,ny,rky,pps,ss)

   !Transform ss back to physical space as u1:
  call spctop_fc(nx,ny,ss,u1,xfactors,yfactors,xtrig,ytrig)
   !Note, u1 here contains d(p1)/dy in physical space

   !(7) Add zonal parts which enforce boundary conditions on zonal mean u:
  do ix=0,nxm1
    do iy=0,ny
      p1(iy,ix)=p1bar(iy)+p1(iy,ix)
      u1(iy,ix)=u1bar(iy)-u1(iy,ix)
    enddo
  enddo

endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !Solve for the mode 2 flow field (p1, u1, v1):

 !(1) remove kd2^2*u2hmean*y from mode 2 PV anomaly to enforce mean zonal
 !    velocity boundary conditions, u2 = u2hmeam at y = ymin & ymax:
dq2dy=kd2sq*u2hmean
do ix=0,nxm1
  do iy=0,ny
    p2(iy,ix)=p2(iy,ix)-dq2dy*yg(iy)
  enddo
enddo

 !(2) FFT p2 and invert to get uncorrected streamfunction p2:
call ptospc_fc(nx,ny,p2,ss,xfactors,yfactors,xtrig,ytrig)
do ky=0,ny
  do kx=0,nxm1
    ss(kx,ky)=qgop2(kx,ky)*ss(kx,ky)
  enddo
enddo

 !(3) Obtain zonal mean streamfunction and zonal velocity:
 !    Extract zonal mean part of p2 from ss and remove it from ss:
sqdn=one/sqrt(dble(nx))
do ky=0,ny
  p2bar(ky)=sqdn*ss(0,ky)
  ss(0,ky)=zero
enddo
 !The sqrt(nx) factor is needed due to a change in the x transform.

 !Calculate a y derivative of the zonal part to get the zonal velocity,
 !apart from the uniform flow u2hmean (note dpdy = 0 on boundaries):
call yderiv_fc(1,ny,rky,p2bar,dpdy)

 !Transform p2bar & dpdy back to physical space (cosine & sine transforms):
call dct(1,ny,p2bar,ytrig,yfactors)
call dst(1,ny, dpdy,ytrig,yfactors)

 !Correct zonal mean streamfunction:
do iy=0,ny
  p2bar(iy)=p2bar(iy)-u2hmean*yg(iy)
enddo

 !Correct zonal mean zonal velocity:
u2bar(0)=u2hmean
do iy=1,nym1
  u2bar(iy)=u2hmean-dpdy(iy)
enddo
u2bar(ny)=u2hmean

 !Transform non-zonal streamfunction back to physical space:
call spctop_fc(nx,ny,ss,p2,xfactors,yfactors,xtrig,ytrig)

 !(4) Do a full transform of p2 at y = ymin and ymax and obtain the
 !    interior field (cppy) that must be subtracted to give p2 = 0
 !    at y = ymin and ymax:
do ix=0,nxm1
  pbot(ix)=p2(0,ix)
  ptop(ix)=p2(ny,ix)
enddo
call forfft(1,nx,pbot,xtrig,xfactors)
call forfft(1,nx,ptop,xtrig,xfactors)

 !Define the (non-zonal, kx > 0) interior semi-spectral field:
do iy=1,nym1
  cppy(iy,0)=zero
enddo
do kx=1,nxm1
  do iy=1,nym1
    cppy(iy,kx)=pbot(kx)*decy2(ny-iy,kx)+ptop(kx)*decy2(iy,kx)
  enddo
enddo
 !Invert using a full transform in x:
call revfft(nym1,nx,cppy,xtrig,xfactors)

 !(5) Remove cppy to obtain the final mode 2 (non-zonal) streamfunction p2:
do ix=0,nxm1
  p2(0, ix)=zero
  p2(ny,ix)=zero
enddo

 !Here we also prepare for the next step:
do ix=0,nxm1
  do iy=1,nym1
    p2(iy,ix)=p2(iy,ix)-cppy(iy,ix)
    v2(iy,ix)=p2(iy,ix)
  enddo
enddo

 !(6) Compute mode 2 velocity field (correct for boundary conditions below):

 !Transform p2 (in v2) to spectral space:
call ptospc_fs(nx,ny,v2,pps,xfactors,yfactors,xtrig,ytrig)

 !Compute d(p2)/dx = ppx spectrally:
call xderiv_fs(nx,ny,hrkx,pps,ppx)

 !Transform ppx back to physical space as v2:
call spctop_fs(nx,ny,ppx,v2,xfactors,yfactors,xtrig,ytrig)

 !Compute d(p2)/dy = ss spectrally:
call yderiv_fs(nx,ny,rky,pps,ss)

 !Transform ss back to physical space as u2:
call spctop_fc(nx,ny,ss,u2,xfactors,yfactors,xtrig,ytrig)
 !Note, u2 here contains d(p2)/dy in physical space

 !(7) Add zonal parts which enforce boundary conditions on zonal mean u:
do ix=0,nxm1
  do iy=0,ny
    p2(iy,ix)=p2bar(iy)+p2(iy,ix)
    u2(iy,ix)=u2bar(iy)-u2(iy,ix)
  enddo
enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !Re-construct layer fields by projecting onto vertical modes:

 !Streamfunction and zonal velocity (add background barotropic flow):
if (tide) then
  ubt=utopo*sin(ftopo*t)
else
  ubt=utopo
endif
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix,1)=vect11*p1(iy,ix)+vect12*p2(iy,ix)
    pp(iy,ix,2)=vect21*p1(iy,ix)+vect22*p2(iy,ix)
    uu(iy,ix,1)=vect11*u1(iy,ix)+vect12*u2(iy,ix)+ubt
    uu(iy,ix,2)=vect21*u1(iy,ix)+vect22*u2(iy,ix)+ubt
  enddo
enddo

 !Meridional velocity (ensure zero edge values):
do ix=0,nxm1
  vv(0,ix,1)=zero
  vv(0,ix,2)=zero
  do iy=1,nym1
    vv(iy,ix,1)=vect11*v1(iy,ix)+vect12*v2(iy,ix)
    vv(iy,ix,2)=vect21*v1(iy,ix)+vect22*v2(iy,ix)
  enddo
  vv(ny,ix,1)=zero
  vv(ny,ix,2)=zero
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

!==================================================

end module     
