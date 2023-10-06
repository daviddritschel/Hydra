program rest
! ---------------------------------------------------------------------|
! |   This routine sets up a state of rest possibly over topography.   |
!----------------------------------------------------------------------|

use constants

implicit none

double precision:: qq(0:ny,0:nxm1),bb(0:ny,0:nxm1)
double precision:: t
integer:: ix,iy

!-----------------------------------------------------------
 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
  open(12,file='topo.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) t,bb
  close(12)
else
  bb=zero
endif

 !Write initial PV distribution to a file 
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)

 !Layer 1 (use beta*y + topographic PV):
do ix=0,nxm1
  do iy=0,ny
    qq(iy,ix)=beta*(ymin+gly*dble(iy))+bb(iy,ix)
  enddo
enddo
write(11,rec=1) zero,qq

 !Layer 2 (use beta*y):
do ix=0,nxm1
  do iy=0,ny
    qq(iy,ix)=beta*(ymin+gly*dble(iy))
  enddo
enddo
write(11,rec=2) zero,qq
close(11)

end program rest
