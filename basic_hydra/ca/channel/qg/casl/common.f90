module common

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use contours
use spectral
use generic

 !Quantities which need to be preserved between recontouring and evolution:

 !Gridded PV fields:
double precision:: qs(0:ny,0:nxm1),qd(0:ny,0:nxm1)

 !For semi-Lagrangian advection: 
double precision:: xig(0:nxm1),yig(0:ny)
double precision,parameter:: xigmax=dble(nx),yigmax=dble(ny) 

 !Initial time of simulation:
double precision:: tinit

 !Time, time step, half time step and time when gridded data are saved:
double precision:: t,dt,hfdt,tgrid,twist

 !Counters for writing direct-access data (gridded and contours):
integer:: igrec,icrec

contains 

!=============================================================
subroutine readcont(iind,iniqs)

! Reads time loop = iind from the contour and residual PV data
! in the cont subdirectory and, optionally, initialises qs if
! iniqs = 1 above.

! Sets the initial time, tinit, from data read from qqsynopsis.asc

implicit none

 !Passed variables:
integer:: iind,iniqs

 !Local variables:
real:: tr4,qdr4(0:ny,0:nxm1)
integer:: iop(nm)
integer:: i,j,ibeg,iend
character(len=3):: pind

write(pind(1:3),'(i3.3)') iind

 !Open contour files:
if (replace) then
   !Here, data from an old run (in restart-*) is used to initialise a new one:
  open(40,file='cont/restart-qqsynopsis.asc',status='old')
  do i=1,iind
    read(40,*) nq,nptq,t,dq,qavg
  enddo
  close(40)

  open(40,file='cont/restart-qqindex'//pind,form='unformatted', &
      & access='direct',status='old',recl=20*nq)
  read(40,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq),iop(1:nq)
  close(40)

  open(40,file='cont/restart-qqnodes'//pind,form='unformatted', &
      & access='direct',status='old',recl=16*nptq)
  read(40,rec=1) xq(1:nptq),yq(1:nptq)
  close(40)

  open(40,file='cont/restart-qqresi'//pind,form='unformatted', &
      & access='direct',status='old',recl=nbytes)
  read(40,rec=1) tr4,qdr4
  qd=dble(qdr4)
  close(40)

else
   !Here, an old run is being continued from a specific time loop, iind:
  open(40,file='cont/qqsynopsis.asc',status='old')
  do i=1,iind
    read(40,*) nq,nptq,t,dq,qavg
  enddo
  close(40)

  open(40,file='cont/qqindex'//pind,form='unformatted', &
      & access='direct',status='old',recl=20*nq)
  read(40,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq),iop(1:nq)
  close(40)

  open(40,file='cont/qqnodes'//pind,form='unformatted', &
      & access='direct',status='old',recl=16*nptq)
  read(40,rec=1) xq(1:nptq),yq(1:nptq)
  close(40)

  open(40,file='cont/qqresi'//pind,form='unformatted', &
      & access='direct',status='old',recl=nbytes)
  read(40,rec=1) tr4,qdr4
  qd=dble(qdr4)
  close(40)
endif

 !Time corresponding to the data read:
tinit=t

!--------------------------------------------------------------------
 !Determine i2q and nextq arrays:
do j=1,nq
  ibeg=i1q(j)
  iend=ibeg+npq(j)-1
  i2q(j)=iend
  do i=ibeg,iend-1
    nextq(i)=i+1
  enddo
  if (iop(j) .eq. 0) then
    nextq(iend)=0
  else
    nextq(iend)=ibeg
  endif
enddo

if (iniqs .eq. 1) then
   !Initialise entire PV field, qs, by combining contours and residual PV:
  call con2grid(qs)
  qs=qs+qd
endif
 
return
end subroutine

!=======================================================================

 !Main end module
end module
