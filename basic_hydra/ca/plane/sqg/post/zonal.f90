program zonal
!  -------------------------------------------------------------------------
!  |   Computes the zonal mean zonal velocity u, <u>, eddy kinetic energy  |
!  |   <(u-<u>)^2+v^2>/2, zonal mean buoyancy, <q>, and eddy enstrophy     |
!  |   <(q-<q>)^2>/2 as a function of y from the grid buoyancy in qq.r4.   |
!  |                                                                       |
!  |   Output is to the unformatted, direct-access file "avg.r4" which     |
!  |   contains single-precision numerical data of the form                |
!  |           t, <u>, <(u-<u>)^2+v^2>/2, <q>, <(q-<q>)^2>/2               |
!  |   with each record corresponding to a given time read from the input  |
!  |   file qq.r4, and with <u>, etc... being real (single-precision)      |
!  |   arrays dimensioned 0:ny (for all grid lines in y, including the     |
!  |   periodic edges 0 & ny, corresponding to y = ymin & ymax).           |
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
double precision:: qa(ny,nx),qq(ny,nx),uu(ny,nx),vv(ny,nx)
 !Spectral fields:
double precision:: qs(nx,ny),pp(nx,ny)

 !Diagnostic quantities:
real:: qqr4(ny,nx)
real:: uavg(0:ny),kavg(0:ny),qavg(0:ny),zavg(0:ny)
real:: t,zfac,hzfac,rmsu,rmsz

integer:: loop,iread,ix,iy

!---------------------------------------------------------------
 !Open file containing buoyancy field:
open(31,file='qq.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

 !Open output files:
open(11,file='eddy.asc',status='replace')
open(22,file='avg.r4',form='unformatted', access='direct', &
                 &  status='replace',recl=4+16*(ny+1))

zfac=1./nx
hzfac=zfac/2.

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) t,qqr4
  if (iread .ne. 0) exit 
  write(*,'(a,f12.5)') ' Processing t = ',t

   !Define the buoyancy anomaly needed for inversion:
  qq=dble(qqr4)
  qa=qq

   !Convert qa to spectral space as qs (note, qa is modified):
  call ptospc(nx,ny,qa,qs,xfactors,yfactors,xtrig,ytrig)

   !Invert buoyancy to get velocity field (uu,vv):
  call main_invert(qs,uu,vv,pp)

   !Obtain zonal averages; first initialise:
  do iy=0,ny-1
    uavg(iy)=0.
    kavg(iy)=0.
    qavg(iy)=0.
    zavg(iy)=0.
  enddo
  
  do ix=1,nx
    do iy=0,ny-1
      uavg(iy)=uavg(iy)+uu(iy+1,ix)
      qavg(iy)=qavg(iy)+qq(iy+1,ix)
    enddo
  enddo

  do iy=0,ny-1
    uavg(iy)=zfac*uavg(iy)
    qavg(iy)=zfac*qavg(iy)
  enddo

  do ix=1,nx
    do iy=0,ny-1
      kavg(iy)=kavg(iy)+(uu(iy+1,ix)-uavg(iy))**2+vv(iy+1,ix)**2
      zavg(iy)=zavg(iy)+(qq(iy+1,ix)-qavg(iy))**2
    enddo
  enddo

  do iy=0,ny-1
    kavg(iy)=hzfac*kavg(iy)
    zavg(iy)=hzfac*zavg(iy)
  enddo
  uavg(ny)=uavg(0)
  kavg(ny)=kavg(0)
  zavg(ny)=zavg(0)

  rmsu=zero
  rmsz=zero
  do iy=0,ny-1
    rmsu=rmsu+kavg(iy)
    rmsz=rmsz+zavg(iy)
  enddo
  rmsu=sqrt(rmsu/dble(ny))
  rmsz=sqrt(rmsz/dble(ny))
  write(11,'(f12.5,2(1x,f16.12))') t,rmsu,rmsz

   !Write diagnostic data for this time:
  write(22,rec=loop) t,uavg,kavg,qavg,zavg
enddo

 !Close output files:
close(11)
close(22)

write(*,*)
write(*,*) ' The results are ready in avg.r4; for information on formatting'
write(*,*) ' see the comments at the top of src/post/zonal.f90'

return
end subroutine

 !End main program
end program
!=======================================================================
