module spectral

use constants
use sta2dfft

 !Maximum x & y wavenumbers:
integer,parameter:: nwx=nx/2,nwy=ny/2

 !Constants to implement spectral forcing:
integer, parameter:: nx2=nx+2,ny2=ny+2
integer, parameter:: nkfmax=nwx*nwy

 !Common arrays, constants:
double precision:: rksq(nx,ny),qdiss(nx,ny)
double precision:: green(nx,ny),flo(nx,ny),fhi(nx,ny)
double precision:: rk(nx,ny)
double precision:: rkx(nx),hrkx(nx)
double precision:: rky(ny),hrky(ny)
double precision:: bety(ny)

 !FFT trig tables:
double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

 !Forcing constants and arrays:
double precision:: dqdt(nx,ny)
double precision:: kxf(nkfmax), kyf(nkfmax), forfac
double precision:: kxfxo(nwx), kyfyo(nwy)
integer:: nkf, nkfxo, nkfyo

double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny))
integer:: kmag(nx,ny),kmax,kmaxred

!=======================================================================!
! From main code: call init_invert               to initialise          !
! then            call main_invert(qq,uu,vv,pp)  to perform inversion   !
!=======================================================================!

contains

!=================================================================
subroutine init_spectral

implicit none

 !Local variables:
double precision:: rkxf(nx),rkyf(ny)
double precision:: scx,rkxmax,scy,rkymax
double precision:: delk,delki,fac,rkmsi,rksqu,snorm
integer:: kx,ky,k,iy

!----------------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define x wavenumbers and filtered x wavenumbers:
rkx(1)=zero
do kx=1,nwx-1
  rkx(kx+1)   =hrkx(2*kx)
  rkx(nx+1-kx)=hrkx(2*kx)
enddo
rkx(nwx+1)=hrkx(nx)

scx=twopi/ellx
rkxmax=scx*dble(nwx)

do kx=1,nx
  rkxf(kx)=rkx(kx)*exp(-36.d0*(rkx(kx)/rkxmax)**36.d0)
  hrkx(kx)=hrkx(kx)*exp(-36.d0*(hrkx(kx)/rkxmax)**36.d0)
enddo

 !Define y wavenumbers and filtered y wavenumbers:
rky(1)=zero
do ky=1,nwy-1
  rky(ky+1)   =hrky(2*ky)
  rky(ny+1-ky)=hrky(2*ky)
enddo
rky(nwy+1)=hrky(ny)

scy=twopi/elly
rkymax=scy*dble(nwy)

do ky=1,ny
  rkyf(ky)=rky(ky)*exp(-36.d0*(rky(ky)/rkymax)**36.d0)
  hrky(ky)=hrky(ky)*exp(-36.d0*(hrky(ky)/rkymax)**36.d0)
enddo

 !Note, these are used when taking derivatives for approximate 
 !de-aliasing (see Hou & Li, J. Nonlinear Sci. 2006).

!-----------------------------------------------------------------------
 !Define squared total wavenumber, diffusion operators, low and hi-pass
 !Butterworth filters (Wireless Engineer, vol. 7, 1930, pp. 536-541),
 !and the filtered squared wavenumber for computing current density:
fac=36.d0/dble(ngridp)
rkmsi=one/max(rkxmax**2,rkymax**2)
do ky=1,ny
  do kx=1,nx
    rksqu=rkx(kx)**2+rky(ky)**2
    rk(kx,ky)=sqrt(rksqu)
    qdiss(kx,ky)=cdamp*(rkmsi*rksqu)**nnu
    flo(kx,ky)=one/(one+(fac*rksqu)**2)
    fhi(kx,ky)=one-flo(kx,ky)
    rksq(kx,ky)=rkxf(kx)**2+rkyf(ky)**2
  enddo
enddo

!-----------------------------------------------------------------------
 !Initialise arrays for computing the spectrum of any field:
delk=sqrt(scx**2+scy**2)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do ky=1,ny
  do kx=1,nx
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

 !Initialise narrow-band spectral forcing if present:
if (forcing) call init_forcing
 
return 
end subroutine

!=================================================================
subroutine init_forcing
! initialises narrow-band spectral forcing

implicit none

double precision:: fran, famp
integer:: k,kx,ky

!------------------------------------
 !Initialise random number generator:
fran=zero
call myrand(-2.d0)
do k=0,iseed
   call myrand(fran)
enddo

!-----------------------------------------------------------------
 !Determine wavenumbers lying in chosen band where forcing occurs:
if (fwid .ge. dble(kfmid)) then
  write(*,*) ' Band extends to 0 wavenumber!  Not allowed!'
  write(*,*) ' *** stopping ***'
  stop
endif

if (fwid .ge. dble(min(nwx-kfmid,nwy-kfmid))) then
  write(*,*) ' Band extends to Nyquist frequency!  Not allowed!'
  write(*,*) ' *** stopping ***'
  stop
endif

 !zero y wavenumber:
nkfxo=0
do kx=1,nx
  if (abs(rkx(kx)-kfmid) .le. fwid+small) then
    nkfxo=nkfxo+1
    kxfxo(nkfxo)=kx
  endif
enddo

 !zero x wavenumber:
nkfyo=0
do ky=1,ny
  if (abs(rky(ky)-kfmid) .le. fwid+small) then
    nkfyo=nkfyo+1
    kyfyo(nkfyo)=ky
  endif
enddo

 !non-zero x and y wavenumbers:
nkf=0
do ky=2,nwy
  do kx=2,nwx
    if (abs(rk(kx,ky)-kfmid) .le. fwid+small) then
      nkf=nkf+1
      kxf(nkf)=kx
      kyf(nkf)=ky
    endif
  enddo
enddo

 !Zero forcing array (only populate non-zero values in subroutine force):
dqdt=zero

 !Forcing amplitude prefactor:
famp=fbase*two**nfamp
forfac=dble(nx*ny/2)*famp/(twopi*fwid)

return
end subroutine

!=================================================================
subroutine nb_forcing(dt)
! fills spectral forcing array dqdt

implicit none

 !Passed variable (time step):
double precision:: dt

 !Local variables:
double precision:: fran1,fran2,ffac,fmag
integer:: k,kxfl,kyfl

!-------------------------------------
 !Modified version of Alvelius (1999):
ffac=forfac*dt

 !Fill array for zero y wavenumber:
do k=1,nkfxo
  kxfl=kxfxo(k)
  call myrand(fran1)
  call myrand(fran2)
  fmag=sqrt(two*ffac*rkx(kxfl)/(fran1**2+fran2**2))
  dqdt(kxfl,1)=fmag*fran1
  dqdt(nx2-kxfl,1)=fmag*fran2
enddo

 !Fill array for zero x wavenumber:
do k=1,nkfyo
  kyfl=kyfyo(k)
  call myrand(fran1)
  call myrand(fran2)
  fmag=sqrt(two*ffac*rky(kyfl)/(fran1**2+fran2**2))
  dqdt(1,kyfl)=fmag*fran1
  dqdt(1,ny2-kyfl)=fmag*fran2
enddo

 !Fill array for non-zero x and y wavenumbers:
do k=1,nkf
  kxfl=kxf(k)
  kyfl=kyf(k)
  call myrand(fran1)
  call myrand(fran2)
  fmag=sqrt(ffac*rk(kxfl,kyfl)/(fran1**2+fran2**2))
  dqdt(kxfl,kyfl)=fmag*fran1
  dqdt(nx2-kxfl,kyfl)=fmag*fran2
  call myrand(fran1)
  call myrand(fran2)
  fmag=sqrt(ffac*rk(kxfl,kyfl)/(fran1**2+fran2**2))
  dqdt(kxfl,ny2-kyfl)=fmag*fran2
  dqdt(nx2-kxfl,ny2-kyfl)=-fmag*fran1
enddo

return
end subroutine

!=================================================================
subroutine myrand(uni)

double precision:: uni
integer:: iseed

iseed=0 !to get rid of unused variable warning
if (abs(uni+two) .lt. small) then
  call srand(0)
else
  uni=two*rand(0)-one
endif

return
end subroutine

!=============================================================
subroutine main_invert(qq,uu,vv,pp)
! Given the PV anomaly qq in spectral space, this routine computes
! the streamfunction pp in spectral space and the velocity field
! (uu,vv) in physical space.

implicit none

 !Passed variables:
double precision:: qq(nx,ny),pp(nx,ny) !Spectral
double precision:: uu(ny,nx),vv(ny,nx) !Physical

 !Local variables:
double precision:: vtmp(nx,ny)
integer:: ix,iy,kx,ky

!--------------------------------
 !Solve for psi (pp):
do ky=1,ny
  do kx=1,nx
    pp(kx,ky)=green(kx,ky)*qq(kx,ky)
  enddo
enddo

 !Get velocity field:
call xderiv(nx,ny,hrkx,pp,vtmp)
call spctop(nx,ny,vtmp,vv,xfactors,yfactors,xtrig,ytrig)

call yderiv(nx,ny,hrky,pp,vtmp)
call spctop(nx,ny,vtmp,uu,xfactors,yfactors,xtrig,ytrig)

 !Copy -uu into uu:
do ix=1,nx
  do iy=1,ny
    uu(iy,ix)=-uu(iy,ix)
  enddo
enddo

return
end subroutine

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
end subroutine

!===================================================================

subroutine spec1d(ss,spec,iopt)
! Computes the 1d spectrum of a spectral field ss and returns the
! result in spec.
! If iopt = 1, the spectral field is multiplied first by K^2
!              and is returned modified!

implicit none

 !Passed variables:
double precision:: ss(nx,ny),spec(0:max(nx,ny))
integer:: iopt

 !Local variables:
integer:: kx,ky,k

!--------------------------------------------------------
if (iopt .eq. 1) then
   !Multiply spectral field by K^2 (filtered, see init_spectral):
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
