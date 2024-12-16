program cross
!  ------------------------------------------------------------------
!  | Computes y = 0 cross sections of scaled fields of b', xi, eta, |
!  | zeta, p & w, as well as z = 0 cross sections of scaled fields  |
!  | of b', xi, eta, zeta, p & delta (= -w_z). To compare with the  |
!  | corresponding unscaled fields in PS3D for example, one should  |
!  | multiply by b' by alpha*f*N, xi by alpha*N, eta by alpha*N,    |
!  | zeta by alpha*f, p by alpha*f^2, w by alpha^2*f/N, and delta   |
!  | by alpha^2, where alpha << 1 is the maximum amplitude of the   |
!  | scaled surface buoyancy, b_0, used in the SQG code. The actual |
!  | surface buoyancy is alpha*f*N*b_0.                             |
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
double precision:: b0(ny,nx),qq(ny,nx)
double precision:: bx(ny,nx),by(ny,nx)
double precision:: px(ny,nx),py(ny,nx)

 !y = 0 cross sections:
double precision:: bby0(0:nz,nx),ppy0(0:nz,nx),wwy0(0:nz,nx)
double precision:: xxy0(0:nz,nx),yyy0(0:nz,nx),zzy0(0:nz,nx)

 !z = 0 cross sections:
double precision:: bbz0(ny,nx),ppz0(ny,nx),ddz0(ny,nx)
double precision:: xxz0(ny,nx),yyz0(ny,nx),zzz0(ny,nx)

 !Spectral fields:
double precision:: bs(nx,ny)
double precision:: wka(nx,ny),wkb(nx,ny)
double precision:: wkf(nx,ny),wkg(nx,ny)
double precision:: wkp(nx,ny),wkm(nx,ny)

 !Work quantities:
double precision:: b0avg,zlin,dz,delt
real:: qqr4(ny,nx)
real:: t

integer:: ny0bytes,loop,iread,iz,ix
character(len=3):: pind

!---------------------------------------------------------
 !Initialise FFTs:
call init_spectral

 !Define unfiltered wavenumber (overwrite rksq):
do ky=1,ny
   do kx=1,nx
      rksq(kx,ky)=sqrt(rkx(kx)**2+rky(ky)**2)
   enddo
enddo
 !Ignore spectral wavenumber corresponding to average:
rksq(1,1)=one

 !Define common spectral arrays:
wkf=filt/(one-exp(-two*rksq*depth))
wkg=wkf/rksq

 !Grid interval in z:
dz=depth/dble(nz)

 !Read average surface buoyancy:
open(10,file='average_qq.asc',status='old')
read(10,*) b0avg
close(10)

!---------------------------------------------------------------
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

 !Open file containing the scaled surface buoyancy field:
open(31,file='bb.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(31,rec=loop) t,qqr4
close(31)
write(*,'(a,f5.1)') ' Processing data at t = ',t
b0=dble(qqr4)
 !Add average buoyancy at surface:
bbz0=b0avg+b0

 !Convert surface bouyancy b0 to spectral space as bs:
call ptospc(nx,ny,b0,bs,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
bs=filt*bs
 !Ensure mean is zero (it is restored below):
bs(1,1)=zero

 !Loop over z and construct (u,v)*grad(b) at each level:
do iz=0,nz
   zlin=dz*dble(iz)
   !Work out b' at height z in spectral space and differentiate:
   !   wka=bs*filt*sinh(rksq*zlin)/sinh(rksq*depth)
   wkp=exp( rksq*(zlin-depth))
   wkm=exp(-rksq*(zlin+depth))
   green=wkf*(wkp-wkm)
   green(1,1)=zero
   wka=bs*green !this is b' in spectral space
   call gradient(wka,bx,by) !this is grad(b) in physical space

   !Save a y = 0 cross section of b' in bby0, as well as cross sections
   !of xi = -b'_x & eta = b'_y in xxy0 & yyy0:
   call spctop(nx,ny,wka,b0,xfactors,yfactors,xtrig,ytrig)
   do ix=1,nx
      bby0(iz,ix)= b0(nwy+1,ix)
      xxy0(iz,ix)=-bx(nwy+1,ix)
      yyy0(iz,ix)=-by(nwy+1,ix)
   enddo

   !Work out phi at height z in spectral space and differentiate:
   !   wka=bs*filt*cosh(rksq*zlin)/(rksq*sinh(rksq*depth))
   green=wkg*(wkp+wkm)
   green(1,1)=zero
   wka=bs*green !this is phi in spectral space
   call gradient(wka,px,py) !this is grad(phi) in physical space
   wkb=-rksq**2*wka !this is zeta in spectral space
   
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

   !Compute (u,v)*grad(b) and filter:
   b0=px*by-py*bx
   call ptospc(nx,ny,b0,wka,xfactors,yfactors,xtrig,ytrig)
   wka=filt*wka
   call spctop(nx,ny,wka,b0,xfactors,yfactors,xtrig,ytrig)
   !Store -(u,v)*grad(b) in a y = 0 cross section:
   do ix=1,nx
      wwy0(iz,ix)=-b0(nwy+1,ix)
   enddo
enddo
 !b0 now contains (u,v)*grad(b) at z = 0.

 !Store xi & eta at z = 0:
xxz0=-bx
yyz0=-by

 !Convert (u,v)*grad(b) at z = 0 to spectral space as bs:
call ptospc(nx,ny,b0,bs,xfactors,yfactors,xtrig,ytrig)

 !Loop over z to add surface contribution:
do iz=0,nz
   zlin=dz*dble(iz)
   green=wkf*(exp(rksq*(zlin-depth))-exp(-rksq*(zlin+depth)))
   green(1,1)=zero
   wka=bs*green
   call spctop(nx,ny,wka,b0,xfactors,yfactors,xtrig,ytrig)
   do ix=1,nx
      wwy0(iz,ix)=wwy0(iz,ix)+b0(nwy+1,ix)
   enddo
enddo

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
write(41,rec=6) t,real(wwy0)
close(41)

 !z = 0 cross sections (excluding w, which is zero):
open(41,file='z0_t'//pind//'.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(41,rec=1) t,real(bbz0)
write(41,rec=2) t,real(xxz0)
write(41,rec=3) t,real(yyz0)
write(41,rec=4) t,real(zzz0)
write(41,rec=5) t,real(ppz0)
write(41,rec=6) t,real(ddz0)
close(41)

write(*,*)
write(*,*) &
     ' y = 0 cross sections of b'', xi, eta, zeta, p & w', &
     ' are in y0_t'//pind//'.r4'
write(*,*) &
     ' z = 0 cross sections of b'', xi, eta, zeta, p & delta', &
     ' are in z0_t'//pind//'.r4'

return
end subroutine diagnose

 !End main program
end program cross
!=======================================================================
