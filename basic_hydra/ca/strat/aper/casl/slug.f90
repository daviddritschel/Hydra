program slug

use constants

! This routine sets up initial buoyancy and vorticity fields for
! casl in a rectangular aperiodic domain.
! The initial condition is a rounded 'slug' buoyancy anomaly, with zero vorticity. 
! Adapted by sek on 5 Sept 2012 from an original f77 code by dgd @ St Andrews

implicit double precision (a-h,o-z)

 !Local vorticity and buoyancy arrays:
double precision:: zz(0:ny,0:nx),bb(0:ny,0:nx)

write(*,'(a,2(f6.2))') ' We consider a domain of height and length: ',elly,ellx
write(*,'(a,2(i6))') ' The y and x grid resolution is: ',ny,nx
gridrat=dble(nx)*elly/(dble(ny)*ellx)
write(*,'(a,2(f6.2))') ' Note, the y:x grid length ratio is: ',gridrat
write(*,*)
write(*,*) ' Enter the width of the higher density zone:'
read(*,*) w

write(*,*) ' Enter the height of the higher density zone:'
read(*,*) h

write(*,*) ' Enter the width of the density transition zone:'
read(*,*) d

 !Set up buoyancy distribution:
do ix=0,nx
  x=glx*dble(ix)
  do iy=0,ny
    y=gly*dble(iy)
    r=sqrt((x/w)**2+(y/h)**2)
    zz(iy,ix)=zero
    bb(iy,ix)=-f12*erfc((r-one)/d)
  enddo
enddo

 !Write vorticity distribution to file:
open(11,file='zz_init.dat',status='unknown')
write(11,'(f3.1)') zero
do ix=0,nx
  do iy=0,ny
    write(11,'(f16.7)') zz(iy,ix)
  enddo
enddo
close(11)

 !Write buoyancy distribution to file:
open(12,file='bb_init.dat',status='unknown')
write(12,'(f3.1)') zero
do ix=0,nx
  do iy=0,ny
    write(12,'(f16.7)') bb(iy,ix)
  enddo
enddo
close(12)

 !Write input data to file:
open(12,file='input_for_slug',status='unknown')
write(12,'(4x,f6.2,a)') w,' ! width of the higher density zone'
write(12,'(4x,f6.2,a)') h,' ! height of the higher density zone'
write(12,'(4x,f6.4,a)') d,' ! width of the density transition zone'
close(12)
      
end program
