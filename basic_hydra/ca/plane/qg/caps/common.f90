module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Define quantities to be preserved between recontouring and evolution:

 !PV residue:
double precision:: qr(ny,nx)

 !PV anomaly evolved in spectral space (note array order):
double precision:: qs(nx,ny)

 !Thermal equilibrium streamfunction (spectral):
double precision,allocatable,dimension(:,:):: ppeq

 !Time stepping parameters:
double precision:: dtmax,tinit,tgrid

 !Parameter for vortex forcing:
double precision:: dnvor

contains 

!=============================================================
subroutine readcont(iind)

! Reads time loop = iind from the contour and residual PV data
! in the cont subdirectory and initialises qs.

! Sets the initial time, tinit, from data read from qqsynopsis.asc

implicit none

 !Passed variable:
integer:: iind

 !Local variables:
double precision:: qa(ny,nx)
real:: tr4,qqr4(ny,nx)
integer:: i,j,ibeg,iend
character(len=3):: pind

!--------------------------------------------------------------------
 !Open contour files to continue from a specific time loop, iind:
open(40,file='cont/qqsynopsis.asc',status='old')
do i=1,iind
  read(40,*) nq,nptq,t,qjump
enddo
close(40)

write(pind(1:3),'(i3.3)') iind

open(40,file='cont/qqindex'//pind,form='unformatted',status='old')
read(40) npq(1:nq),i1q(1:nq),indq(1:nq)
close(40)

open(40,file='cont/qqnodes'//pind,form='unformatted',status='old')
read(40) xq(1:nptq),yq(1:nptq)
close(40)

open(40,file='cont/qqresi.r4',form='unformatted',status='old', &
                          & access='direct',recl=nbytes)
read(40,rec=iind) tr4,qqr4
close(40)
qr=dble(qqr4)

 !Time corresponding to the data read:
tinit=t

!--------------------------------------------------------------------
 !Set contour ending indices (i2q) and nextq array:
do j=1,nq
  ibeg=i1q(j)
  iend=ibeg+npq(j)-1
  i2q(j)=iend
  do i=ibeg,iend-1
    nextq(i)=i+1
  enddo
  nextq(iend)=ibeg
enddo

!--------------------------------------------------------------------
 !Initialise entire PV field, qs, by combining contours and residual PV:

 !Convert PV contours to gridded values (qa):
call con2grid(qa)

 !Add residual in qr:
qa=qa+qr

 !Convert qa to spectral space as qs (note, qa is modified):
call ptospc(nx,ny,qa,qs,xfactors,yfactors,xtrig,ytrig)
 
return
end subroutine

!=======================================================================

end module
