module spectral

use constants
use sta2dfft

 !Maximum x & y wavenumber:
integer,parameter:: nwx=nx/2,nwxm1=nwx-1,nwxp1=nwx+1
integer,parameter:: nwy=ny/2,nwym1=nwy-1,nwyp1=nwy+1

 !Common arrays, constants:
double precision:: green(nx,ny),rksq(nx,ny),diff(nx,ny),rksqf(nx,ny)
double precision:: rkx(nx),hrkx(nx)
double precision:: rky(ny),hrky(ny)
double precision:: bety(ny)

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(nx,ny),kmax,kmaxred

!====================================================================!
! From main code: call init_invert            to initialise          !
! then            call main_invert(zz,uu,vv)  to perform inversion   !
!====================================================================!

contains

!===========================
subroutine init_spectral

implicit none

 !Local variables:
double precision:: rkxf(nx),rkyf(ny)
double precision:: scx,rkxmax,scy,rkymax
double precision:: delk,delki,snorm
integer:: kx,ky,k,iy

!----------------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define x wavenumbers:
rkx(1)=zero
do kx=1,nwxm1
  rkx(kx+1)   =hrkx(2*kx)
  rkx(nx+1-kx)=hrkx(2*kx)
enddo
rkx(nwx+1)=hrkx(nx)
scx=twopi/ellx
rkxmax=scx*dble(nwx)

 !Define y wavenumbers:
rky(1)=zero
do ky=1,nwym1
  rky(ky+1)   =hrky(2*ky)
  rky(ny+1-ky)=hrky(2*ky)
enddo
rky(nwy+1)=hrky(ny)
scy=twopi/elly
rkymax=scy*dble(nwy)

 !Define squared total wavenumber:
do ky=1,ny
  do kx=1,nx
    rksq(kx,ky)=rkx(kx)**2+rky(ky)**2
  enddo
enddo

 !Define approximate de-aliasing filter (ref: Hou & Li, J. Nonlinear Sci. 
 !2006) and apply in taking derivatives:
do kx=1,nx
  rkxf(kx)=rkx(kx)*exp(-36.d0*(rkx(kx)/rkxmax)**36.d0)
  hrkx(kx)=hrkx(kx)*exp(-36.d0*(hrkx(kx)/rkxmax)**36.d0)
enddo
do ky=1,ny
  rkyf(ky)=rky(ky)*exp(-36.d0*(rky(ky)/rkymax)**36.d0)
  hrky(ky)=hrky(ky)*exp(-36.d0*(hrky(ky)/rkymax)**36.d0)
enddo

 !Define squared total filtered wavenumber used for computing j:
do ky=1,ny
  do kx=1,nx
    rksqf(kx,ky)=rkxf(kx)**2+rkyf(ky)**2
  enddo
enddo

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

 !Define beta*y:
do iy=1,ny
  bety(iy)=beta*(gly*dble(iy-1)+ymin)
enddo
 
return 
end subroutine

!=========================================
subroutine main_invert(qq,uu,vv,pp,ppflag)

implicit none

 !Passed variables:
double precision:: qq(ny,nx),uu(ny,nx),vv(ny,nx),pp(ny,nx)
logical:: ppflag

 !Local variables:
double precision:: ss(nx,ny),vtmp(nx,ny)
integer:: ix,iy,kx,ky

!--------------------------------
 !Solve for psi (pp):
do ix=1,nx
  do iy=1,ny
    pp(iy,ix)=qq(iy,ix)
  enddo
enddo

call ptospc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)

do ky=1,ny
  do kx=1,nx
    ss(kx,ky)=green(kx,ky)*ss(kx,ky)
  enddo
enddo

 !Get velocity field:
call xderiv(nx,ny,hrkx,ss,vtmp)
call spctop(nx,ny,vtmp,vv,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,ss,vtmp)
call spctop(nx,ny,vtmp,uu,xfactors,yfactors,xtrig,ytrig)

 !Copy -uu into uu:
do ix=1,nx
  do iy=1,ny
    uu(iy,ix)=-uu(iy,ix)
  enddo
enddo

 !Get psi in physical space if required:
if (ppflag) call spctop(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!=================================================================
subroutine lorentz(apot,qqsrc)
! Computes the magnetic PV forcing, qqsrc, defined by
!      S_b = B_0*dj/dx - J(A,j)
! where j = -Lap{A} is the current density.  This routine also
! computes dA/dx & dA/dy and returns them in the arrays aax & aay.
! Note: apot = A

implicit none

 !Passed arrays:
double precision:: apot(ny,nx),qqsrc(ny,nx)

 !Local variables:
double precision:: aax(ny,nx),aay(ny,nx)
double precision:: jjx(nx,ny),jjy(nx,ny)
double precision:: ss(nx,ny),vtmp(nx,ny)
integer:: ix,iy,kx,ky

!---------------------------------------------------
 !Get A in spectral space (use jjx as a work array):
do ix=1,nx
  do iy=1,ny
    jjx(iy,ix)=apot(iy,ix)
  enddo
enddo

call ptospc(nx,ny,jjx,ss,xfactors,yfactors,xtrig,ytrig)

 !Get derivatives of A (in temporary array ss):
call xderiv(nx,ny,hrkx,ss,vtmp)
call spctop(nx,ny,vtmp,aax,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,ss,vtmp)
call spctop(nx,ny,vtmp,aay,xfactors,yfactors,xtrig,ytrig)

 !Define j in spectral space (use ss again):
do ky=1,ny
  do kx=1,nx
    ss(kx,ky)=rksqf(kx,ky)*ss(kx,ky)
  enddo
enddo

 !Get derivatives of j:
call xderiv(nx,ny,hrkx,ss,vtmp)
call spctop(nx,ny,vtmp,jjx,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,ss,vtmp)
call spctop(nx,ny,vtmp,jjy,xfactors,yfactors,xtrig,ytrig)

 !Compute source term, S_b:
do ix=1,nx
  do iy=1,ny
    qqsrc(iy,ix)=(b0+aay(iy,ix))*jjx(iy,ix)-aax(iy,ix)*jjy(iy,ix)
  enddo
enddo

return
end subroutine

!===============================================================
subroutine magfield(apot,aax,aay)
! Computes the magnetic field anomaly in (aay,-aax) where 
! aax = dA/dx & aay = dA/dy.
! Note: apot = A

implicit none

 !Passed arrays:
double precision:: apot(ny,nx),aax(ny,nx),aay(ny,nx)

 !Local variables:
double precision:: ss(nx,ny),vtmp(nx,ny)
integer:: ix,iy

!---------------------------------------------------
 !Get A in spectral space (use vtmp as a work array):
do ix=1,nx
  do iy=1,ny
    vtmp(iy,ix)=apot(iy,ix)
  enddo
enddo

call ptospc(nx,ny,vtmp,ss,xfactors,yfactors,xtrig,ytrig)

 !Get derivatives of A:
call xderiv(nx,ny,hrkx,ss,vtmp)
call spctop(nx,ny,vtmp,aax,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,ss,vtmp)
call spctop(nx,ny,vtmp,aay,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!===============================================================
subroutine current(apot,jj)
! Computes the current density j from the magnetic potential A
! Note: apot = A

implicit none

 !Passed arrays:
double precision:: apot(ny,nx),jj(ny,nx)

 !Local variables:
double precision:: ss(nx,ny)
integer:: ix,iy,kx,ky

!----------------------------------------------------
 !Get A in spectral space (use vtmp as a work array):
do ix=1,nx
  do iy=1,ny
    jj(iy,ix)=apot(iy,ix)
  enddo
enddo

call ptospc(nx,ny,jj,ss,xfactors,yfactors,xtrig,ytrig)

 !Define j in spectral space (use ss again):
do ky=1,ny
  do kx=1,nx
    ss(kx,ky)=rksqf(kx,ky)*ss(kx,ky)
  enddo
enddo

 !FFT back to physical space:
call spctop(nx,ny,ss,jj,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!===================================================================

subroutine diffuse(apot,aad,hfdt)
! Diffuses a field A and defines a diffusion operator diff needed in
! subroutine spread below.  hfdt = dt/2 where dt is the time step.
! Note: apot = A

implicit none

 !Passed variables:
double precision:: apot(ny,nx),aad(ny,nx)
double precision:: hfdt

 !Local variables:
double precision:: ss(nx,ny)
double precision:: alpha
integer:: ix,iy,kx,ky

!--------------------------------
 !Get A in spectral space:
do ix=1,nx
  do iy=1,ny
    aad(iy,ix)=apot(iy,ix)
  enddo
enddo

call ptospc(nx,ny,aad,ss,xfactors,yfactors,xtrig,ytrig)

 !Apply diffusion and define diff operator:
alpha=hfdt*eta
do ky=1,ny
  do kx=1,nx
    diff(kx,ky)=one/(one+alpha*rksq(kx,ky))
    ss(kx,ky)=ss(kx,ky)*diff(kx,ky)*(one-alpha*rksq(kx,ky))
  enddo
enddo

 !Return diffused A to physical space:
call spctop(nx,ny,ss,aad,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!===================================================================

subroutine spread(ff)
! Spreads a field F using the diffusion operator diff defined above
! in subroutine diffuse.

implicit none

 !Passed array:
double precision:: ff(ny,nx)

 !Local variables:
double precision:: ss(nx,ny)
integer:: kx,ky

!-----------------------------------------------------
 !Get F in spectral space:
call ptospc(nx,ny,ff,ss,xfactors,yfactors,xtrig,ytrig)

 !Apply diff operator:
do ky=1,ny
  do kx=1,nx
    ss(kx,ky)=ss(kx,ky)*diff(kx,ky)
  enddo
enddo

 !Return diffused F to physical space:
call spctop(nx,ny,ss,ff,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!===================================================================

subroutine spec1d(rvar,spec,iopt)
! Computes the 1d spectrum of a field rvar which is
! periodic in x and y, and returns the result in spec.
! If iopt = 1, the spectral field is multiplied first by K^2.
!   *** Warning: rvar is modified by this routine.***

implicit none

 !Passed variables:
double precision:: rvar(ny,nx),spec(0:max(nx,ny))
integer:: iopt

 !Local variables:
double precision:: ss(nx,ny)
integer:: kx,ky,k

!--------------------------------------------------------
 !Transform rvar to spectral space:
call ptospc(nx,ny,rvar,ss,xfactors,yfactors,xtrig,ytrig)

if (iopt .eq. 1) then
   !Multiply transformed field by K^2:
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=ss(kx,ky)*rksq(kx,ky)
    enddo
  enddo
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
end subroutine

!===================================================================

end module     
