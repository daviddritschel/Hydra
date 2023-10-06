program gauss
! Initialises a Gaussian vortex with zero divergence in gradient-wind
! (cyclo-geostrophic) balance.

 ! Import constants and parameters:
use constants
 ! Import spectral module:
use spectral

implicit none

double precision,parameter:: glu=twopi/dble(ngu)
double precision:: zod0(ngu/2),zod1(ngu/2),zod2(ngu/2)
double precision:: zev0(0:ngu/2),zev1(0:ngu/2),zev2(0:ngu/2)
double precision:: za(ngu,ngu),zz(ng,ng)
double precision:: uu(ng,ng),vv(ng,ng),gg(ng,ng)
double precision:: hh(ng,ng),qq(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng)
double precision:: eps,a,zmax,aisq,x,y,zavg,hfac
double precision:: bx0,bz0
integer:: ix,iy,ngh,nguf,miy,mix,mixm1

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

write(*,*) ' We take zeta = z_max*exp(-0.5*(r/a)^2), where'
write(*,*) ' z_max = 1.270747*z_0 and z_0 = eps*f.'
write(*,*)
write(*,*) ' Enter the Rossby number, eps:'
read(*,*) eps
write(*,*) ' Enter the radius, a:'
read(*,*) a
write(*,*) ' We take the horizontal magnetic field B = (Bx,0) with'
write(*,*) ' Bx = Bx_0*H/h.  Enter Bx_0 (constant):'
read(*,*) bx0
write(*,*) ' We take Bz = B_z0.  Enter Bz_0 (constant):'
read(*,*) bz0

!Choose zmax so that velocity at r = R is equal to z0*R/2 = u0:
zmax=eps*cof/(2.d0*(1.d0-exp(-0.5d0)))
aisq=one/(two*a**2)

do ix=1,ngu
  x=glu*dble(ix-1)-pi
  do iy=1,ngu
    y=glu*dble(iy-1)-pi
    za(iy,ix)=zmax*exp(-aisq*(x**2+y**2))
  enddo
enddo

!------------------------------------------------------------------------
 !Average the vorticity field in za to the coarser grid (ng,ng):
ngh=ngu
do while (ngh .gt. ng)
  nguf=ngh
  ngh=ngh/2
   !Perform nine-point averaging:
  do iy=1,ngh
    miy=2*iy
    zod2(iy)=za(miy-1,nguf)
    zev2(iy)=za(miy,nguf)
  enddo
  zev2(0)=za(nguf,nguf)
  do ix=1,ngh
    mix=2*ix
    mixm1=mix-1
    do iy=1,ngh
      miy=2*iy
      zod1(iy)=za(miy-1,mixm1)
      zod0(iy)=za(miy-1,mix)
      zev1(iy)=za(miy,mixm1)
      zev0(iy)=za(miy,mix)
    enddo
    zev1(0)=zev1(ngh)
    zev0(0)=zev0(ngh)
    do iy=1,ngh
      za(iy,ix)=0.0625d0*(zev0(iy)+zev0(iy-1)+zev2(iy)+zev2(iy-1)) &
              & +0.125d0*(zev1(iy)+zev1(iy-1)+zod0(iy)+zod2(iy)) &
              &   +0.25d0*zod1(iy)
    enddo
    do iy=1,ngh
      zod2(iy)=zod0(iy)
      zev2(iy)=zev0(iy)
    enddo
    zev2(0)=zev0(0)
  enddo
enddo

 !Calculate and remove average zz:
zavg=zero
do ix=1,ng
  do iy=1,ng
    zavg=zavg+za(iy,ix)
  enddo
enddo
zavg=zavg/dble(ng*ng)

do ix=1,ng
  do iy=1,ng
    zz(iy,ix)=za(iy,ix)-zavg
  enddo
enddo

!----------------------------------------------------------------
! Invert vorticity for streamfunction:
call ptospc(ng,ng,zz,wkb,xfactors,yfactors,xtrig,ytrig)
wka=rlap*wkb !wka = psi in spectral space
wkb=filt*wkb !spectrally-truncated vorticity
call spctop(ng,ng,wkb,zz,xfactors,yfactors,xtrig,ytrig)
! Get non-diverent velocity field:
call gradient(wka,vv,uu)
! uu has the wrong sign, but this is convenient for the next step:

! Compute gamma_l:
call jacob(uu,vv,gg)
gg=two*gg
! Convert to spectral space as wkb:
call ptospc(ng,ng,gg,wkb,xfactors,yfactors,xtrig,ytrig)

! Obtain dimensionless height anomaly:
hfac=one/csq
wka=hfac*(cof*wka-rlap*wkb)
call spctop(ng,ng,wka,hh,xfactors,yfactors,xtrig,ytrig)

! Filter gamma_l and return to physical space:
wkb=filt*wkb
call spctop(ng,ng,wkb,gg,xfactors,yfactors,xtrig,ytrig)

! Compute linearised PV:
qq=zz-cof*hh

!----------------------------------------------------------------
 !Write linearised PV:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

 !Write height field (purely diagnostic):
open(11,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,hh
close(11)

qq=bx0/(one+hh)

 !Write zero velocity divergence:
hh=zero
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,hh
close(11)

 !Write acceleration divergence:
open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,gg
close(11)

 !Write uniform magnetic field components:
open(11,file='bx_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

hh=zero
open(11,file='by_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,hh
close(11)

hh=bz0
open(11,file='bz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,hh
close(11)

end program
