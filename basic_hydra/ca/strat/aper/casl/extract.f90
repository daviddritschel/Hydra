program extract

!  ------------------------------------------------------------------
!  |   Extracts bouyancy and vorticity from data in bb.dat          |
!  |   and zz.dat at a given frame, writing them to bb_ex.dat and   |
!  |   zz_ex.dat.                                                   |
!  ------------------------------------------------------------------

use constants

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: bb(0:ny,0:nx),zz(0:ny,0:nx)

!----------------------------------------------------------------

write(*,*)
write(*,*) ' Enter the frame to extract: '
read(*,*) kfr

 !Open input data file:
open(15,file='zz.dat',status='old')
open(16,file='bb.dat',status='old')

if (kfr .gt. 1) then
  do i=1,(kfr-1)*(nxp1*nyp1+1)
    read(15,*)
    read(16,*) 
  enddo
endif

 !Read field and process:
read(15,*) t
read(16,*) t
write(*,'(a,f12.5)') 'Extracting t= ',t
do ix=0,nx
  do iy=0,ny
    read(15,*) zz(iy,ix)
    read(16,*) bb(iy,ix)
  enddo
enddo
close(15)
close(16)

 !Write vorticity distribution to file:
open(11,file='zz_ex.dat',status='unknown')
write(11,'(f3.1)') zero
do ix=0,nx
  do iy=0,ny
    write(11,'(f16.7)') zz(iy,ix)
  enddo
enddo
close(11)

 !Write buoyancy distribution to file:
open(12,file='bb_ex.dat',status='unknown')
write(12,'(f3.1)') zero
do ix=0,nx
  do iy=0,ny
    write(12,'(f16.7)') bb(iy,ix)
  enddo
enddo
close(12)

!=========================================================

end program
