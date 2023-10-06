module spectral

use constants
use sta2dfft

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Maximum x & y wavenumber:
integer,parameter:: nwx=nx/2,nwxm1=nwx-1,nwxp1=nwx+1
integer,parameter:: nwy=ny/2,nwym1=nwy-1,nwyp1=nwy+1

 !Common arrays, constants:
double precision:: green(nx,ny),rksq(nx,ny),diff(nx,ny),deal(nx,ny)
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

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

!----------------------

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
scy=pi/elly
rkymax=scy*dble(ny)

 !Define squared total wavenumber and dealising filter:
f49=4.d0/9.d0
do ky=1,ny
  do kx=1,nx
    rksq(kx,ky)=rkx(kx)**2+rky(ky)**2
    if ((rkx(kx)/rkxmax)**2+(rky(ky)/rkymax)**2 .lt. f49) then
      deal(kx,ky)=one
    else
      deal(kx,ky)=zero
    endif
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

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(ny,nx),uu(ny,nx),vv(ny,nx),pp(ny,nx)
logical:: ppflag
 !Local arrays:
double precision:: ss(nx,ny),vtmp(nx,ny)

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
subroutine lorentz(aa,qqsrc)
! Computes the magnetic PV forcing, qqsrc, defined by
!      S_b = B_0*dj/dx - J(A,j)
! where j = -Lap{A} is the current density.  This routine also
! computes dA/dx & dA/dy and returns them in the arrays aax & aay.

! ==> de-aliasing is used here - and only here - to compute S_b <==

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: aa(ny,nx),qqsrc(ny,nx)
 !Local arrays:
double precision:: aax(ny,nx),aay(ny,nx)
double precision:: jjx(nx,ny),jjy(nx,ny)
double precision:: ss(nx,ny),vtmp(nx,ny)

!---------------------------------------------------
 !Get A in spectral space (use jjx as a work array):
do ix=1,nx
  do iy=1,ny
    jjx(iy,ix)=aa(iy,ix)
  enddo
enddo

call ptospc(nx,ny,jjx,ss,xfactors,yfactors,xtrig,ytrig)

 !Apply de-aliasing (see initialisation above):
!do ky=1,ny
!  do kx=1,nx
!    ss(kx,ky)=deal(kx,ky)*ss(kx,ky)
!  enddo
!enddo

 !Get derivatives of A (in temporary array ss):
call xderiv(nx,ny,hrkx,ss,vtmp)
call spctop(nx,ny,vtmp,aax,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,ss,vtmp)
call spctop(nx,ny,vtmp,aay,xfactors,yfactors,xtrig,ytrig)

 !Define j in spectral space (use ss again):
do ky=1,ny
  do kx=1,nx
    ss(kx,ky)=rksq(kx,ky)*ss(kx,ky)
  enddo
enddo

 !Get derivatives of j:
call xderiv(nx,ny,hrkx,ss,vtmp)
call spctop(nx,ny,vtmp,jjx,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,ss,vtmp)
call spctop(nx,ny,vtmp,jjy,xfactors,yfactors,xtrig,ytrig)

 !Compute ***fully de-aliased*** source term, S_b:
do ix=1,nx
  do iy=1,ny
    qqsrc(iy,ix)=(b0+aay(iy,ix))*jjx(iy,ix)-aax(iy,ix)*jjy(iy,ix)
  enddo
enddo

return
end subroutine

!===============================================================
subroutine magfield(aa,aax,aay)
! Computes the magnetic field anomaly in (aay,-aax) where 
! aax = dA/dx & aay = dA/dy

! No de-aliasing is done here as the output is used only diagnostically

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: aa(ny,nx),aax(ny,nx),aay(ny,nx)
 !Local arrays:
double precision:: ss(nx,ny),vtmp(nx,ny)

!---------------------------------------------------
 !Get A in spectral space (use vtmp as a work array):
do ix=1,nx
  do iy=1,ny
    vtmp(iy,ix)=aa(iy,ix)
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

!===================================================================

subroutine diffuse(aa,aad,hfdt)
! Diffuses a field A and defines a diffusion operator diff needed in
! subroutine spread below.  hfdt = dt/2 where dt is the time step.

 !Passed arrays:
double precision:: aa(ny,nx),aad(ny,nx)

 !Local array:
double precision:: ss(nx,ny)

!--------------------------------
 !Get A in spectral space:
do ix=1,nx
  do iy=1,ny
    aad(iy,ix)=aa(iy,ix)
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

 !Passed arrays:
double precision:: ff(ny,nx)

 !Local array:
double precision:: ss(nx,ny)

!--------------------------------
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

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: rvar(ny,nx),spec(0:max(nx,ny))
 !Local array:
double precision:: wks(nx,ny)

 !Transform rvar to spectral space:
call ptospc(nx,ny,rvar,wks,xfactors,yfactors,xtrig,ytrig)

if (iopt .eq. 1) then
   !Multiply transformed field by K^2:
  do ky=1,ny
    do kx=1,nx
      wks(kx,ky)=wks(kx,ky)*rksq(kx,ky)
    enddo
  enddo
endif

do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(1,1)
spec(k)=spec(k)+f14*wks(1,1)**2

 !y-independent mode:
do kx=2,nx
  k=kmag(kx,1)
  spec(k)=spec(k)+f12*wks(kx,1)**2
enddo

 !x-independent mode:
do ky=2,ny
  k=kmag(1,ky)
  spec(k)=spec(k)+f12*wks(1,ky)**2
enddo

 !All other modes:
do ky=2,ny
  do kx=2,nx
    k=kmag(kx,ky)
    spec(k)=spec(k)+wks(kx,ky)**2
  enddo
enddo

return
end subroutine

!===================================================================

end module     
