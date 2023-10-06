program tracercont
! ======================================================================= !
! Extracts the tracer contours at a specified time and writes them to the !
! file pvcont.dat, which can be plotted using ptc.py in the scripts       !
! subdirectory.                                                           !
!                                                                         !
! Written 3 November 2017 by D G Dritschel @ St Andrews                   !
! ======================================================================= !

use contours

implicit none

 !Declarations:
double precision:: t
integer:: iind,i,j,ibeg
character(len=3):: pind

!-----------------------------------------------------------------
write(*,*) ' Initialising contours module...'
call init_contours

write(*,*)
write(*,*) ' Time frame to extract tracer contours?'
read(*,*) iind
write(pind(1:3),'(i3.3)') iind

 !Read data:
open(40,file='contours/qqsynopsis.asc',status='old')
do i=1,iind
  read(40,*) nq,nptq,t,qjump
enddo
close(40)

write(*,'(a,f12.5)') ' This corresponds to t = ',t

 !Read contour indices:
open(40,file='contours/qqindex'//pind,form='unformatted', &
      access='direct',status='old',recl=12*nq)
read(40,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq)
close(40)

 !Read contour nodes:
open(40,file='contours/qqnodes'//pind,form='unformatted', &
      access='direct',status='old',recl=16*nptq)
read(40,rec=1) xq(1:nptq),yq(1:nptq)
close(40)

 !Work out number of active contours and nodes, na & npta:
j=nq
do while (indq(j) .eq. 9999)
  j=j-1
enddo
na=j
npta=i1q(j)+npq(j)-1

 !Open output file and write data:
open(55,file='contours/pvcont.dat',status='replace')
write(55,'(i6,1x,i8,1x,f12.5)') nq-na,nptq-npta,t
ibeg=1
do j=na+1,nq
  write(55,'(i6,1x,i8,1x,i1)') npq(j),ibeg,1
  ibeg=ibeg+npq(j)
enddo
do i=npta+1,nptq
  write(55,'(2(1x,f14.11))') xq(i),yq(i)
enddo
close(55)

end program tracercont
