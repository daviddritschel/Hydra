program dambreak

use constants

! This routine sets up initial buoyancy and vorticity fields for
! casl in a rectangular aperiodic domain.
! The setup is to produce a rectangular dam-break with a smoothed discontiuity
! in the buoyancy field. The vorticity can be set to contain random noise.

! Adapted by sek on 5 Sept 2012 from original F77 code by dgd @ St Andrews
implicit double precision (a-h,o-z)

 !Local vorticity and buoyancy arrays:
double precision:: zz(0:ny,0:nx),bb(0:ny,0:nx)

 !Read in parameter values for setup:     
write(*,'(a,2(f6.2))') ' We consider a domain of height and length: ',elly,ellx
write(*,'(a,2(i6))') ' The y and x grid resolution is: ',ny,nx
gridrat=dble(nx)*elly/(dble(ny)*ellx)
write(*,'(a,2(f6.2))') ' Note, the y:x grid length ratio is: ',gridrat
write(*,*)
write(*,*) ' Enter the width of the higher density zone:'
read(*,*) w

write(*,*) ' Enter the height of the higher density zone:'
read(*,*) h              !Not used if b1 is commented out below.

write(*,*) ' Enter the width of the density transition zone:'
read(*,*) d
 
write(*,*) ' Max. amplitude of vorticity noise:'
read(*,*) amp

write(*,*) ' Random seed for noise:'
read(*,*) seed
iseed=int(seed)

 !Initialize random # generator:
do i=1,iseed
  uni=rand(0)
enddo

 !Set up buoyancy distribution:
do ix=0,nx
  x=glx*dble(ix)
  b0=-f12*erfc((x-w)/d)
  do iy=0,ny
    y=gly*dble(iy)
    b1=f12*erfc((y-h)/d)
    zz(iy,ix)=amp*(two*rand(0)-one)
    bb(iy,ix)=b0!*b1     !Include b1 if vertical structure is required 
  enddo
enddo

 !Write vorticity distribution to file:
open(20,file='zz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

 !Write buoyancy distribution to file:
open(20,file='bb_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,bb
close(20)

 !Write input data to file:
!open(12,file='input_for_dam-break',status='unknown')
!write(12,'(4x,f6.2,a)') w,' ! width of the higher density zone'
!write(12,'(4x,f6.2,a)') h,' ! height of the higher density zone'
!write(12,'(4x,f6.4,a)') d,' ! width of the density transition zone'
!write(12,'(1x,f9.5,a)') amp,' ! Max. amplitude of vorticity noise'
!write(12,'(3x,i7,a)') iseed,' ! Random seed for vorticity noise'
!close(12)

      
contains 

!==========================================================================
double precision function rand(i)
 !Returns a random number: i is any integer

call random_number(r)
rand=r

return
end function

!==========================================================================
end program


