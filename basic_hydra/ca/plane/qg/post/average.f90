program average
!  -------------------------------------------------------------------------
!  |      Computes the time and zonal mean zonal velocity u and PV q       |
!  -------------------------------------------------------------------------

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use spectral

implicit none

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(ny,nx),uu(ny,nx),vv(ny,nx),wk(ny,nx)
double precision:: qs(nx,ny),ps(nx,ny)

 !Diagnostic quantities:
real:: qqr4(ny,nx)
real:: uavg(0:ny),qavg(0:ny),uzon(0:ny),qzon(0:ny)
real:: t1,t2,t,zfac,hzfac,yg

integer:: loop,iread,ix,iy,navg

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open file containing PV field:
open(31,file='qq.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

zfac=1./nx
hzfac=zfac/2.

write(*,*) ' Time period to average, t_1 & t_2?'
read(*,*) t1,t2
t1=t1-1.e-6
t2=t2+1.e-6

!---------------------------------------------------------------
 !Read data and process:
loop=0
navg=0
uavg=0.
qavg=0.
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) t,qqr4
  if (iread .ne. 0) exit 
  if ((t2-t)*(t-t1) .gt. 0.) then
    write(*,'(a,f12.5)') ' Processing t = ',t

     !Define the PV anomaly needed for inversion:
    do ix=1,nx
      do iy=1,ny
        qq(iy,ix)=dble(qqr4(iy,ix))-bety(iy)
      enddo
    enddo

     !Convert to spectral space for inversion:
    wk=qq
    call ptospc(nx,ny,wk,qs,xfactors,yfactors,xtrig,ytrig)

     !Invert spectral PV (qs) to get velocity field (uu,vv):
    call main_invert(qs,uu,vv,ps)

     !Obtain zonal averages; first initialise:
    do iy=0,ny-1
      uzon(iy)=0.
      qzon(iy)=0.
    enddo
  
    do ix=1,nx
      do iy=0,ny-1
        uzon(iy)=uzon(iy)+uu(iy+1,ix)
        qzon(iy)=qzon(iy)+qq(iy+1,ix)
      enddo
    enddo

    do iy=0,ny-1
      uzon(iy)=zfac*uzon(iy)
      qzon(iy)=zfac*qzon(iy)
    enddo
    uzon(ny)=uzon(0)
    qzon(ny)=qzon(0)

    navg=navg+1
    uavg=uavg+uzon
    qavg=qavg+qzon
  endif
enddo

 !Complete time average:
zfac=1.0/float(navg)
uavg=zfac*uavg
qavg=zfac*qavg

 !Add back beta*y:
do iy=0,ny-1
  qavg(iy)=qavg(iy)+bety(iy+1)
enddo
qavg(ny)=qavg(0)+beta*elly

 !Write data:

 !Open output file:
open(21,file='uavg.asc',status='replace')
open(22,file='qavg.asc',status='replace')
do iy=0,ny
  yg=gly*float(iy)-hly
  write(21,'(2(1x,f15.10))') uavg(iy),yg
  write(22,'(2(1x,f15.10))') qavg(iy),yg
enddo
close(21)
close(22)

write(*,*)
write(*,*) ' The results are ready in uavg.asc & qavg.asc.  Plot with'
write(*,*)
write(*,*) ' plotcol uavg.asc'
write(*,*) ' plotcol qavg.asc'
write(*,*)

return
end subroutine

 !End main program
end program average
!=======================================================================
