program rest

! Sets up a rest state consisting of PV varying as
!                2*Omega*z

use constants 

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: slat(ng),hh(ng,nt),qq(ng,nt)

!------------------------------------------------------------
 !Define sin(latitude):
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  slat(j)=sin(rlat)
enddo

 !Define hh (used for writing h & d) and PV:
do i=1,nt
  do j=1,ng
    hh(j,i)=zero
    qq(j,i)=fpole*slat(j)
  enddo
enddo

 !Write initial height field:
open(20,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

 !Write initial divergence field:
open(20,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

 !Write thermal equilibrium height field:
open(20,file='hequil.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)


end program
