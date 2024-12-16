module spectral

use constants
use sta2dfft

 !Maximum x & y wavenumbers:
integer,parameter:: nwx=nx/2,nwy=ny/2

 !Common arrays, constants:
double precision:: rksq(nx,ny),qdiss(nx,ny)
double precision:: green(nx,ny),vorop(nx,ny)
double precision:: filt(nx,ny),flo(nx,ny),fhi(nx,ny)
double precision:: rkx(nx),hrkx(nx)
double precision:: rky(ny),hrky(ny)

 !For diffusing a tracer (optional):
double precision,allocatable,dimension(:,:):: tdiss

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(nx,ny),kmax,kmaxred

!=======================================================================!
! From main code: call init_invert               to initialise          !
! then            call main_invert(qq,uu,vv,pp)  to perform inversion   !
!=======================================================================!

contains

!===========================
subroutine init_spectral

implicit none

 !Local variables:
double precision:: scx,rkxmax,scy,rkymax,rkfsq
double precision:: delk,delki,fac,rkmsi,rksqu,snorm
integer:: kx,ky,k,iy

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
snorm=-two*depth
do ky=1,ny
  do kx=1,nx
    fac=exp(snorm*sqrt(rkx(kx)**2+rky(ky)**2))
    green(kx,ky)=(one+fac)/(one-fac)
  enddo
enddo

green(1,1)=zero
do kx=2,nx
  green(kx,1)=green(kx,1)/rkx(kx)
enddo
do ky=2,ny
  green(1,ky)=green(1,ky)/rky(ky)
enddo
do ky=2,ny
  do kx=2,nx
    green(kx,ky)=green(kx,ky)/sqrt(rkx(kx)**2+rky(ky)**2)
  enddo
enddo

!-----------------------------------------------------------------------
 !Define squared total wavenumber, diffusion operators, low and hi-pass
 !Butterworth filters (Wireless Engineer, vol. 7, 1930, pp. 536-541), etc:
fac=36.d0/dble(ngridp)
rkfsq=dble(ngridp)/9.d0
rkmsi=one/max(rkxmax**2,rkymax**2)
do ky=1,ny
  do kx=1,nx
    rksqu=rkx(kx)**2+rky(ky)**2
    qdiss(kx,ky)=cdamp*(rkmsi*rksqu)**nnu
    if (rksqu .gt. rkfsq) then
      filt(kx,ky)=zero
      flo(kx,ky)=zero
      fhi(kx,ky)=zero
      rksq(kx,ky)=zero
      green(kx,ky)=zero
      vorop(kx,ky)=zero
    else
      filt(kx,ky)=one
      flo(kx,ky)=one/(one+(fac*rksqu)**2)
      fhi(kx,ky)=one-flo(kx,ky)
      rksq(kx,ky)=rksqu
      vorop(kx,ky)=-green(kx,ky)*rksqu
    endif
  enddo
enddo

if (tracer) then
   !Define tracer diffusion operator:
  allocate(tdiss(nx,ny))
  do ky=1,ny
    do kx=1,nx
      tdiss(kx,ky)=kappa*(rkx(kx)**2+rky(ky)**2)
    enddo
  enddo
endif

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
    k=nint(sqrt(rksq(kx,ky))*delki)
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

return 
end subroutine init_spectral

!=========================================
subroutine main_invert(qq,uu,vv,pp)
! Given the buoyancy anomaly qq in spectral space, this routine computes
! the streamfunction pp in spectral space and the velocity field
! (uu,vv) in physical space.

implicit none

 !Passed variables:
double precision:: qq(nx,ny),pp(nx,ny) !Spectral
double precision:: uu(ny,nx),vv(ny,nx) !Physical

 !Local variable:
double precision:: vtmp(nx,ny) !Spectral

!--------------------------------
 !Solve for psi (pp):
pp=green*qq

 !Get velocity field:
call xderiv(nx,ny,hrkx,pp,vtmp)
call spctop(nx,ny,vtmp,vv,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,pp,vtmp)
call spctop(nx,ny,vtmp,uu,xfactors,yfactors,xtrig,ytrig)

 !Switch sign of uu:
uu=-uu

return
end subroutine main_invert

!=================================================================
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

subroutine spec1d(ss,spec)
! Computes the 1d spectrum of a spectral field ss and returns the
! result in spec.

implicit none

 !Passed variables:
double precision:: ss(nx,ny),spec(0:max(nx,ny))

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

end module spectral
