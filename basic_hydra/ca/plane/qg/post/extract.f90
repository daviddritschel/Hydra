program extract
!=========================================================================!
!  Computes the streamfunction at a specified time and writes the data
!  for imaging with dataview.

!  Written 26/7/2018 by dgd @ WC
!=========================================================================!

 !Import constants and parameters:
use parameters
use constants
 !Import St Andrews FFTs:
use sta2dfft

implicit none

 !Maximum x & y wavenumbers:
integer,parameter:: nwx=nx/2,nwy=ny/2

 !Arrays needed for inversion:
double precision:: green(nx,ny),rkx(nx),rky(ny),bety(ny)
double precision:: ss(nx,ny)

 !Physical fields:
double precision:: qq(ny,nx),pp(ny,nx)
real:: qqr4(ny,nx),ppr4(ny,nx)

 !FFT trig tables:
double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

 !Other local quantities:
double precision:: hrkx(nx),hrky(ny)
real:: t
integer:: kx,ky,ix,iy,k1
character(len=4):: pind

!--------------------------------------------------------------
write(*,*) ' Enter the time to extract q & psi:'
read(*,*) t
k1=nint(t/tgsave)+1
write(pind,'(i4.4)') k1

!----------------------------------------------------------------------
 !Initialise FFTs for inversion:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

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
  rky(ny+1-ky)=hrky(2*ky)
enddo
rky(nwy+1)=hrky(ny)

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

!---------------------------------------------------------------
 !Open file containing PV field:
open(31,file='qq.r4',form='unformatted',access='direct', &
     status='old',recl=nbytes)

 !Open output files:
open(21,file='qq'//pind//'.r4',form='unformatted',access='direct', &
     status='replace',recl=nbytes)
open(22,file='pp'//pind//'.r4',form='unformatted',access='direct', &
     status='replace',recl=nbytes)

!---------------------------------------------------------------
 !Read data and process:
read(31,rec=k1) t,qqr4
write(*,'(a,f12.5)') ' Processing t = ',t

if (beffect) then
   !Define the PV anomaly needed for inversion:
  do ix=1,nx
    qq(:,ix)=dble(qqr4(:,ix))-bety
  enddo
else
  qq=dble(qqr4)   
endif
 
 !Convert qq to spectral space (as pp temporarily):
call ptospc(nx,ny,qq,ss,xfactors,yfactors,xtrig,ytrig)

 !Invert PV to get the streamfunction:
ss=green*ss
call spctop(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)
ppr4=real(pp)

 !Write data for this time:
write(21,rec=1) t,qqr4
write(22,rec=1) t,ppr4

 !Close all files:
close(21)
close(22)
close(31)

write(*,*)
write(*,*) ' To view the results, type'
write(*,'(a,i5,1x,i5,a)') ' dataview qq'//pind//'.r4 -ndim ',nx,ny,' &'
write(*,'(a,i5,1x,i5,a)') ' dataview pp'//pind//'.r4 -ndim ',nx,ny,' &'

 !End main program
end program
!=======================================================================




