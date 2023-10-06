program zonal
!  -------------------------------------------------------------------------
!  |   Computes the zonal mean zonal velocity u, <u>, eddy kinetic energy  |
!  |   <(u-<u>)^2+v^2>/2, zonal mean PV, <q>, and eddy enstrophy           |
!  |   <(q-<q>)^2>/2 as a function of y from the grid PV in qq.r4.         |
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
double precision:: qq(ny,nx),uu(ny,nx),vv(ny,nx),pp(ny,nx)

 !Diagnostic quantities:
real:: qqr4(ny,nx)
real:: uavg(0:ny),kavg(0:ny),qavg(0:ny),zavg(0:ny)
real:: t,zfac,hzfac

integer:: loop,iread,ix,iy

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open file containing PV field:
open(31,file='qq.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

 !Open output file:
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

   !Define the PV anomaly needed for inversion:
  do ix=1,nx
    do iy=1,ny
      pp(iy,ix)=dble(qqr4(iy,ix))-bety(iy)
    enddo
  enddo

   !Convert to spectral space for inversion:
  call ptospc(nx,ny,pp,qq,xfactors,yfactors,xtrig,ytrig)

   !Invert PV to get velocity field (uu,vv):
  call main_invert(qq,uu,vv,pp)

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

   !Add back beta*y:
  do iy=0,ny-1
    qavg(iy)=qavg(iy)+bety(iy+1)
  enddo
  qavg(ny)=qavg(0)+beta*elly

   !Write diagnostic data for this time:
  write(22,rec=loop) t,uavg,kavg,qavg,zavg
enddo

 !Close output file:
close(22)

write(*,*)
write(*,*) ' The results are ready in avg.r4; for information on formatting'
write(*,*) ' see the comments at the top of src/post/zonal.f90'

return
end subroutine

 !End main program
end program
!=======================================================================
