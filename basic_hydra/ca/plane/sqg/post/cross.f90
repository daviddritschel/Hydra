program cross
!  ------------------------------------------------------------------
!  | Computes y = 0 cross sections of scaled fields of b', xi, eta, |
!  | zeta, p, u, v & w, and z = 0 cross sections of scaled fields   |
!  | of b', xi, eta, zeta, p, u, v & u_x + v_y. To compare with the |
!  | corresponding unscaled fields in PS3D for example, one should  |
!  | multiply by b' by alpha*f*N, xi by alpha*N, eta by alpha*N,    |
!  | zeta by alpha*f, p by alpha*f^2, u & v by alpha*f, w by        |
!  | alpha^2*f^2/N, and delta = u_x + v_y by alpha^2*f, where       |
!  | alpha << 1 is the maximum amplitude of the scaled surface      |
!  | buoyancy, b_0, used in the SQG code. The actual surface        |
!  | buoyancy is alpha*f*N*b_0.                                     |
!  |                                                                |  
!  | Outputs y0_tnnn.r4, cross sections at y = 0 at time "nnn", as  |  
!  | well as z0_tnnn.r4, cross sections at z = 0 at time "nnn".     |
!  |                                                                |  
!  | Completed 11 April 2024 by D G Dritschel @ St Andrews          |  
!  ------------------------------------------------------------------

 !Import constants, parameters and arrays needed for FFTs:
use spectral

implicit none

double precision:: tp
integer:: nz,kx,ky

!---------------------------------------------------------
write(*,*) ' Enter time to process:'
read(*,*) tp

 !Work out number of z grid intervals for an isotropic grid in
 !coordinates x & Nz/f:
nz=nint(dble(nx)*depth/ellx)
write(*,'(a,i4,a)') ' Using ',nz,' vertical intervals.'
write(*,*)

 !Read data and process:
call diagnose(tp,nz)

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose(tp,nz)

implicit none

 !Passed variables:
double precision:: tp
integer:: nz

 !Physical fields:
double precision:: b0(ny,nx),z0(ny,nx)
double precision:: bx(ny,nx),by(ny,nx)
double precision:: px(ny,nx),py(ny,nx)
double precision:: zx(ny,nx),zy(ny,nx)

 !y = 0 cross sections:
double precision:: bby0(0:nz,nx),ppy0(0:nz,nx),wwy0(0:nz,nx)
double precision:: xxy0(0:nz,nx),yyy0(0:nz,nx),zzy0(0:nz,nx)
double precision:: uuy0(0:nz,nx),vvy0(0:nz,nx)

 !z = 0 cross sections:
double precision:: bbz0(ny,nx),ppz0(ny,nx),ddz0(ny,nx)
double precision:: xxz0(ny,nx),yyz0(ny,nx),zzz0(ny,nx)
double precision:: uuz0(ny,nx),vvz0(ny,nx)

 !Spectral fields:
double precision:: rk(nx,ny),bs(nx,ny)
double precision:: wka(nx,ny),wkb(nx,ny)
double precision:: wkf(nx,ny),wkg(nx,ny)
double precision:: wkp(nx,ny),wkm(nx,ny)

 !Work quantities:
double precision:: zlin(0:nz),bavg(0:nz)
double precision:: dz,b0avg,delt
real:: qqr4(ny,nx)
real:: t

integer:: ny0bytes,loop,iread,iz,ix
character(len=3):: pind

!-----------------------------------------------------------------
 !Initialise FFTs:
call init_spectral

 !Define unfiltered squared wavenumber K^2 (rksq):
do ky=1,ny
   do kx=1,nx
      rksq(kx,ky)=rkx(kx)**2+rky(ky)**2
   enddo
enddo
rksq(1,1)=zero

 !Define unfiltered wavenumber K (rk):
rk=sqrt(rksq)
 !Ignore spectral wavenumber corresponding to average:
rk(1,1)=one

 !Define common spectral arrays:
wkf=filt/(one-exp(-two*rk*depth))
wkg=wkf/rk

 !Grid interval in z:
dz=depth/dble(nz)

 !Read average surface buoyancy:
open(10,file='average_qq.asc',status='old')
read(10,*) b0avg
close(10)

 !Define average buoyancy at any height z:
do iz=0,nz
   zlin(iz)=dz*dble(iz) !this is z + D (ranges from 0 to D = NH/f)
   bavg(iz)=b0avg*zlin(iz)/depth
enddo
 !depth stands for D in the notes

!-----------------------------------------------------------------
 !Find appropriate record number to read for time chosen:
delt=0.1d0*tgsave
if (tp < delt) then
   loop=1
else
   !Read ene-ens.asc to find record number:
   open(21,file='ene-ens.asc',status='old')
   t=zero
   loop=0
   do while (t+delt < tp)
      read(21,*) t
      loop=loop+1
   enddo
   close(21)
endif

!-----------------------------------------------------------------
 !Open file containing the scaled surface buoyancy field:
open(31,file='bb.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(31,rec=loop) t,qqr4
close(31)
write(*,'(a,f5.1)') ' Processing data at t = ',t
b0=dble(qqr4)
 !Add average buoyancy at surface:
bbz0=b0avg+b0
 !b0 below does not include this average; add as needed below.

 !Convert surface bouyancy b0 to spectral space as bs:
call ptospc(nx,ny,b0,bs,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
bs=filt*bs
 !Ensure mean is zero for consistency with calculations below:
bs(1,1)=zero

!-----------------------------------------------------------------
 !Loop over z and construct (u,v)*grad(b) & (u,v)*grad(zeta) at each z:
do iz=0,nz
   !Work out b' at height z in spectral space and differentiate:
   !   wka=bs*sinh(K*(z+D))/sinh(K*D) will store b'
   wkp=exp( rk*(zlin(iz)-depth))
   wkm=exp(-rk*(zlin(iz)+depth))
   green=wkf*(wkp-wkm)      !this is sinh(K*(z+D))/sinh(K*D)
   wka=bs*green             !this is b' in spectral space
   call gradient(wka,bx,by) !this is grad(b') in physical space

   !Save a y = 0 cross section of b' in bby0, as well as cross sections
   !of xi = -b'_x & eta = -b'_y in xxy0 & yyy0:
   call spctop(nx,ny,wka,b0,xfactors,yfactors,xtrig,ytrig)
   do ix=1,nx
      bby0(iz,ix)= bavg(iz)+b0(nwy+1,ix) !Need to add average b here
      xxy0(iz,ix)=-bx(nwy+1,ix)
      yyy0(iz,ix)=-by(nwy+1,ix)
   enddo

   !Work out phi at height z in spectral space and differentiate:
   !   wka=bs*cosh(K*(z+D))/(K*sinh(K*D)) will store phi
   green=wkg*(wkp+wkm)      !this is cosh(K*(z+D))/(K*sinh(K*D))
   wka=bs*green             !this is phi in spectral space
   call gradient(wka,px,py) !this is grad(phi) in physical space
   wkb=-rksq*wka            !this is zeta in spectral space
   call gradient(wkb,zx,zy) !this is grad(zeta) in physical space
   
   !Save a y = 0 cross section of phi in ppy0:
   call spctop(nx,ny,wka,ppz0,xfactors,yfactors,xtrig,ytrig)
   do ix=1,nx
      ppy0(iz,ix)=ppz0(nwy+1,ix)
   enddo

   !Save a y = 0 cross section of zeta in zzy0:
   call spctop(nx,ny,wkb,zzz0,xfactors,yfactors,xtrig,ytrig)
   do ix=1,nx
      zzy0(iz,ix)=zzz0(nwy+1,ix)
   enddo

   !Save y = 0 cross section of u & v in uuy0 & vvy0:
   do ix=1,nx
      uuy0(iz,ix)=-py(nwy+1,ix)
      vvy0(iz,ix)= px(nwy+1,ix)
   enddo

   !Compute (u,v)*grad(zeta) and filter (store in z0):
   z0=px*zy-py*zx
   call ptospc(nx,ny,z0,wka,xfactors,yfactors,xtrig,ytrig)
   wka=filt*wka
   call spctop(nx,ny,wka,z0,xfactors,yfactors,xtrig,ytrig)

   !Compute (u,v)*grad(b') and filter (store in b0):
   b0=px*by-py*bx
   call ptospc(nx,ny,b0,wka,xfactors,yfactors,xtrig,ytrig)
   wka=filt*wka
   call spctop(nx,ny,wka,b0,xfactors,yfactors,xtrig,ytrig)

   !Store (u,v)*grad(b') in a y = 0 cross section:
   do ix=1,nx
      wwy0(iz,ix)=b0(nwy+1,ix)
   enddo
enddo
 !b0 & z0 now contain (u,v)*grad(b') & (u,v)*grad(zeta) at z = 0.

 !Store xi & eta at z = 0:
xxz0=-bx
yyz0=-by

 !Store u & v at z = 0:
uuz0=-py
vvz0= px

 !Convert (u,v)*grad(b') at z = 0 to spectral space as bs:
call ptospc(nx,ny,b0,bs,xfactors,yfactors,xtrig,ytrig)
bs(1,1)=zero

 !Loop over z to add surface contribution to w throughout domain:
do iz=0,nz
   wkp=exp( rk*(zlin(iz)-depth))
   wkm=exp(-rk*(zlin(iz)+depth))
   green=wkf*(wkp-wkm)      !this is sinh(K*(z+D))/sinh(K*D)
   wka=bs*green             !this is -b'_t in spectral space
   call spctop(nx,ny,wka,b0,xfactors,yfactors,xtrig,ytrig)
   !Store a y = 0 cross section of w:
   do ix=1,nx
      wwy0(iz,ix)=b0(nwy+1,ix)-wwy0(iz,ix)
   enddo
enddo

 !Complete calculation of the surface divergence delta:
green=rksq*wkg*(wkp+wkm) !this is K*cosh(K*D)/sinh(K*D)
wka=bs*green             !this is zeta_t in spectral space
call spctop(nx,ny,wka,b0,xfactors,yfactors,xtrig,ytrig)
 !Store divergence at z = 0:
ddz0=-z0-b0

!-----------------------------------------------------------------
 !Open output files and write data:
ny0bytes=4*(nx*(nz+1)+1)
write(pind,'(i3.3)') nint(tp)

 !y = 0 cross sections:
open(41,file='y0_t'//pind//'.r4',form='unformatted', &
      access='direct',status='replace',recl=ny0bytes)
write(41,rec=1) t,real(bby0)
write(41,rec=2) t,real(xxy0)
write(41,rec=3) t,real(yyy0)
write(41,rec=4) t,real(zzy0)
write(41,rec=5) t,real(ppy0)
write(41,rec=6) t,real(uuy0)
write(41,rec=7) t,real(vvy0)
write(41,rec=8) t,real(wwy0)
close(41)

 !z = 0 cross sections (excluding w, which is zero):
open(41,file='z0_t'//pind//'.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(41,rec=1) t,real(bbz0)
write(41,rec=2) t,real(xxz0)
write(41,rec=3) t,real(yyz0)
write(41,rec=4) t,real(zzz0)
write(41,rec=5) t,real(ppz0)
write(41,rec=6) t,real(uuz0)
write(41,rec=7) t,real(vvz0)
write(41,rec=8) t,real(ddz0)
close(41)

!-----------------------------------------------------------------
write(*,*)
write(*,*) &
     ' y = 0 cross sections of b'', xi, eta, zeta, p, u, v & w', &
     ' are in y0_t'//pind//'.r4'
write(*,*) &
     ' z = 0 cross sections of b'', xi, eta, zeta, p, u, v & delta', &
     ' are in z0_t'//pind//'.r4'

return
end subroutine diagnose

 !End main program
end program cross
!=======================================================================
