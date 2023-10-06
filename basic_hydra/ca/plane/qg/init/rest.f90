program rest
! Initialises a flow at rest

use constants

implicit none

double precision:: qq(ny,nx)
integer:: ix,iy

if (beffect) then
  do ix=1,nx
    do iy=1,ny
      qq(iy,ix)=beta*(ymin+gly*dble(iy-1))
    enddo
  enddo
else
  qq=zero
endif

 !Write initial PV (beta*y):
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

 !Relax back to rest:
qq=zero
open(11,file='psieq.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program
