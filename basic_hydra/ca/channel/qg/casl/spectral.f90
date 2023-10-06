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
double precision:: decy(nym1,nxm1)
double precision:: green(0:nxm1,0:ny)
double precision:: rkx(0:nxm1),hrkx(nx),rky(ny)

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

double precision:: frkx(0:nxm1),frky(ny)
double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(0:nxm1,0:ny),kmax,kmaxred

!==========================================================================!
! From main code: call init_invert                   to initialise         !
! then            call main_invert(qq,uu,vv,pp,zz)   to perform inversion  !
!==========================================================================!

contains

!===========================
subroutine init_spectral

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

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
 !Define inversion operator (green):
green(0,0)=-one/kdsq
do kx=1,nxm1
  rksq=rkx(kx)**2
  green(kx,0)=-one/(rksq+kdsq)
enddo
do ky=1,ny
  rksq=rky(ky)**2
  green(0,ky)=-one/(rksq+kdsq)
enddo
do ky=1,ny
  do kx=1,nxm1
    rksq=rkx(kx)**2+rky(ky)**2
    green(kx,ky)=-one/(rksq+kdsq)
  enddo
enddo

!----------------------------------------------------------------------
 !Hyperbolic function used to correct boundary conditions in inversion:
do kx=1,nxm1
  fac=sqrt(rkx(kx)**2+kdsq)*elly
  div=one/(one-exp(-two*fac))
  do iy=1,nym1
    argm=fac*(one-yh1(iy))
    argp=fac*(one+yh1(iy))
    decy(iy,kx)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

return
end subroutine

!========================================
subroutine main_invert(qq,uu,vv,pp,zz)

 !Input:  PV field (qq)
 !Output: Velocity field (uu,vv), streamfunction (pp), and
 !        relative vorticity (zz).

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(0:ny,0:nxm1),pp(0:ny,0:nxm1)
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: zz(0:ny,0:nxm1)
 !Local arrays:
double precision:: cppy(nym1,0:nxm1)
double precision:: vi(ny,0:nxm1)
double precision:: ss(0:nxm1,0:ny),pps(0:nxm1,ny),ppx(0:nxm1,ny)
double precision:: pbot(0:nxm1),ptop(0:nxm1)
double precision:: ppbar(0:ny),uubar(0:ny)
double precision:: dpdy(ny)

!---------------------------------------------------------------------
 !Store PV anomaly in pp temporarily:
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix)=qq(iy,ix)-bety(iy)
  enddo
enddo
 !Also prepare zz to contain relative vorticity at end of routine:
zz=pp

 !FFT pp -> ss and invert to get uncorrected streamfunction:
call ptospc_fc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)
ss=green*ss

 !Extract zonal mean part of pp from ss and remove it from ss:
sqdn=one/sqrt(dble(nx))
do ky=0,ny
  ppbar(ky)=sqdn*ss(0,ky)
  ss(0,ky)=zero
enddo
 !The sqrt(nx) factor is needed due to a change in the x transform.

 !Calculate a y derivative of the zonal part to get the zonal mean
 !zonal velocity, -dpdy (note dpdy = 0 on boundaries):
call yderiv_fc(1,ny,rky,ppbar,dpdy)

 !Transform ppbar & uubar back to physical space (cosine & sine transforms):
call dct(1,ny,ppbar,ytrig,yfactors)
call dst(1,ny, dpdy,ytrig,yfactors)

 !Correct zonal mean zonal velocity and ensure it has mean ubar:
const=ubar+dnyi*sum(dpdy(1:nym1))
uubar(0)=const
do iy=1,nym1
  uubar(iy)=const-dpdy(iy)
enddo
uubar(ny)=const

 !Transform non-zonal streamfunction ss back to physical space as pp:
call spctop_fc(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

 !Do a full transform of pp at y = ymin and ymax and obtain the
 !interior field (cppy) that must be subtracted to give pp = 0
 !at y = ymin and ymax:
do ix=0,nxm1
  pbot(ix)=pp(0,ix)
  ptop(ix)=pp(ny,ix)
enddo
call forfft(1,nx,pbot,xtrig,xfactors)
call forfft(1,nx,ptop,xtrig,xfactors)

 !Define the (non-zonal, kx > 0) interior semi-spectral field:
do iy=1,nym1
  cppy(iy,0)=zero
enddo
do kx=1,nxm1
  do iy=1,nym1
    cppy(iy,kx)=pbot(kx)*decy(ny-iy,kx)+ptop(kx)*decy(iy,kx)
  enddo
enddo
 !Invert using a full transform in x:
call revfft(nym1,nx,cppy,xtrig,xfactors)

 !Remove cppy to obtain the final non-zonal streamfunction pp:
do ix=0,nxm1
  pp(0, ix)=zero
  pp(ny,ix)=zero
enddo

 !Here we also prepare to calculate v = dpsi/dx (here vi):
do ix=0,nxm1
  do iy=1,nym1
    pp(iy,ix)=pp(iy,ix)-cppy(iy,ix)
    vi(iy,ix)=pp(iy,ix)
  enddo
enddo

 !Transform non-zonal streamfunction to spectral space:
call ptospc_fs(nx,ny,vi,pps,xfactors,yfactors,xtrig,ytrig)

 !Compute dpsi/dx = ppx spectrally:
call xderiv_fs(nx,ny,hrkx,pps,ppx)

 !Transform ppx back to physical space as vi (interior part of v):
call spctop_fs(nx,ny,ppx,vi,xfactors,yfactors,xtrig,ytrig)

 !Compute dpsi/dy = ss spectrally:
call yderiv_fs(nx,ny,rky,pps,ss)

 !Transform ss back to physical space as uu:
call spctop_fc(nx,ny,ss,uu,xfactors,yfactors,xtrig,ytrig)
 !Note, uu here contains dpsi/dy in physical space

 !Add zonal parts from above:
do ix=0,nxm1
  do iy=0,ny
    pp(iy,ix)=ppbar(iy)+pp(iy,ix)
    uu(iy,ix)=uubar(iy)-uu(iy,ix)
  enddo
enddo

 !Meridional velocity, now including zero edge values:
do ix=0,nxm1
  vv(0,ix)=zero
  do iy=1,nym1
    vv(iy,ix)=vi(iy,ix)
  enddo
  vv(ny,ix)=zero
enddo

 !Compute relative vorticity:
zz=zz+kdsq*pp

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
