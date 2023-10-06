program random_topo
! |-----------------------------------------------------------|
! |      This routine sets up random bottom topography.       |
! |      Optionally, one can include a linear slope in y.     |
! |-----------------------------------------------------------|

use constants
use sta2dfft

implicit none

 !Doubled domain in y:
integer,parameter:: nye=2*ny
 !Maximum x & y wavenumbers in domain doubled in y (and periodic):
integer,parameter:: nwx=nx/2,nwy=nye/2

 !Common arrays, constants:
double precision:: bb(nye,nx)
double precision:: ss(nx,nye)
double precision:: br(0:ny,nx)
double precision:: hrkx(nx),hrky(nye),rkx(nx),rky(nye)
double precision:: xtrig(2*nx),ytrig(2*nye)
double precision:: fac,q,s,phix,phiy,rms,p,rk0
double precision:: amp,cx,sx,cy,sy,slope

integer:: xfactors(5),yfactors(5)
integer:: kx,ky,kxc,kyc,k,i,ix,iy
integer, dimension(:), allocatable :: seed

!-----------------------------------------------------------------------
 !Initialise random number generator:
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=iseed
do i=1,iseed
  call random_seed(put=seed)
enddo

!----------------------------------------------------------------------
 !Initialise FFTs:
call init2dfft(nx,nye,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define x wavenumbers:
rkx(1)=zero
do kx=1,nwx-1
  rkx(kx+1)   =hrkx(2*kx)
  rkx(nx+1-kx)=hrkx(2*kx)
enddo
rkx(nwx+1)=hrkx(nx)

 !Define y wavenumbers:
rky(1)=zero
do ky=1,nwy-1
  rky(ky+1)   =hrky(2*ky)
  rky(nye+1-ky)=hrky(2*ky)
enddo
rky(nwy+1)=hrky(nye)

!-----------------------------------------------------
 !Compute random field with a given variance spectrum:
write(*,*) ' We consider a topography spectrum of the form'
write(*,*) ' c k^{2p+1} exp[-p(k/k_0)^2]. Enter p > 0 and k_0:'
read(*,*) p,rk0
write(*,*) ' Enter the rms value of the field (determines c):'
read(*,*) rms
write(*,*) ' This field is superposed on top of a linear sloping'
write(*,*) ' topography in y with zero mean height. Enter the slope:'
read(*,*) slope

fac=one/rk0**2
q=p/two
do ky=1,nwy+1
  do kx=1,nwx+1
    s=fac*(rkx(kx)**2+rky(ky)**2)
    ss(kx,ky)=s**q*exp(-q*s)
  enddo
enddo

! Apply to generate full spectrum:
do ky=2,nwy
  kyc=nye+2-ky
  do kx=2,nwx
    kxc=nx+2-kx
    call random_number(s)
    phix=twopi*s-pi
    call random_number(s)
    phiy=twopi*s-pi
    cx=cos(phix)
    sx=sin(phix)
    cy=cos(phiy)
    sy=sin(phiy)
    amp=ss(kx,ky)
    ss(kx ,ky )=amp*cx*cy
    ss(kxc,ky )=amp*sx*cy
    ss(kx, kyc)=amp*cx*sy
    ss(kxc,kyc)=amp*sx*sy
  enddo
enddo

ky=1
do kx=2,nwx
  kxc=nx+2-kx
  call random_number(s)
  phix=twopi*s-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cx
  ss(kxc,ky )=amp*sx
enddo

kx=1
do ky=2,nwy
  kyc=nye+2-ky
  call random_number(s)
  phiy=twopi*s-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cy
  ss(kx, kyc)=amp*sy
enddo

ky=nwy+1
do kx=2,nwx
  kxc=nx+2-kx
  call random_number(s)
  phix=twopi*s-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cx
  ss(kxc,ky )=amp*sx
enddo

kx=nwx+1
do ky=2,nwy
  kyc=nye+2-ky
  call random_number(s)
  phiy=twopi*s-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=ss(kx,ky)
  ss(kx ,ky )=amp*cy
  ss(kx, kyc)=amp*sy
enddo

! Remove mean and Nyquist frequency:
ss(1,1)=zero
ss(nwx+1,nwy+1)=zero

! Transform to physical space:
call spctop(nx,nye,ss,bb,xfactors,yfactors,xtrig,ytrig)

! Work out rms:
fac=sqrt(sum(bb**2)/dble(nx*nye))

! Renormalise field:
fac=rms/fac
bb=fac*bb

! Add uniform sloping part in y:
do ix=1,nx
  do iy=0,ny
    br(iy,ix)=bb(iy+1,ix)+slope*(gly*dble(iy)-hly)
  enddo
enddo

! Write data to a file:
open(11,file='topo.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,br
close(11)

end program random_topo
