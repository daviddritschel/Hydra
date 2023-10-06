program psiq
!=========================================================================!
!  Computes the streamfunction on a selected coarse grid over a selected
!  time range.

!  Written 26/7/2018 by dgd @ WC
!=========================================================================!

 !Import constants and parameters:
use parameters
use constants
 !Import St Andrews FFTs:
use sta2dfft

implicit none

 !Inversion grid/coarse grid ratio:
integer,parameter:: mgc=16,nxc=nx/mgc,nyc=ny/mgc

 !Maximum x & y wavenumbers:
integer,parameter:: nwx=nx/2,nwy=ny/2

 !Arrays needed for inversion:
double precision:: green(nx,ny),rkx(nx),rky(ny),bety(ny)
double precision:: ss(nx,ny)

 !Physical fields:
double precision:: qq(ny,nx),pp(ny,nx)
real:: qqr4(ny,nx),qqc(nyc,nxc),ppc(nyc,nxc)

 !FFT trig tables:
double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

 !Other local quantities:
double precision:: hrkx(nx),hrky(ny)
double precision:: t1,t2,delt
real:: t
integer:: kx,ky,ix,iy,ixc,iyc
integer:: k1,k2,km,k,loop

!--------------------------------------------------------------
write(*,*) ' Enter the range of times to process, t_1 and t_2:'
read(*,*) t1,t2
write(*,*) ' Enter the time interval between processed frames, delta_t:'
read(*,*) delt
k1=nint(t1/tgsave)+1
k2=nint(t2/tgsave)+1
km=nint(delt/tgsave)

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

 !Open output file:
open(22,file='pq.r4',form='unformatted',access='direct', &
     status='replace',recl=4+8*nyc*nxc)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do k=k1,k2,km
  read(31,rec=k) t,qqr4
  write(*,'(a,f12.5)') ' Processing t = ',t

  if (beffect) then
     !Define the PV anomaly needed for inversion:
    do ix=1,nx
      qq(:,ix)=dble(qqr4(:,ix))-bety
    enddo
  else
    qq=dble(qqr4)   
  endif

   !Convert qq to spectral space (as ss temporarily):
  call ptospc(nx,ny,qq,ss,xfactors,yfactors,xtrig,ytrig)

   !Invert PV to get the streamfunction:
  ss=green*ss
  call spctop(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

   !Sample data and write to output file:
  do ixc=1,nxc
    ix=1+mgc*(ixc-1)
    do iyc=1,nyc
      iy=1+mgc*(iyc-1)
      ppc(iyc,ixc)=real(pp(iy,ix))
      qqc(iyc,ixc)=qqr4(iy,ix)
    enddo
  enddo

   !Write diagnostic data for this time:
  loop=loop+1
  write(22,rec=loop) t,ppc,qqc
enddo

 !Close output file:
close(22)

write(*,*)
write(*,*) ' The results are ready in pq.r4.'

 !End main program
end program
!=======================================================================




