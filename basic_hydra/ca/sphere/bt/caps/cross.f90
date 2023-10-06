program cross

! Program which extracts a cross section of constant latitude or 
! longitude from data previously created by genfg

implicit none

integer:: ntx,nty
character(len=8):: infile

!-----------------------------------------------------------------
write(*,*) ' Enter the file in the fine subdirectory to use:'
read(*,*) infile

open(44,file='fine/'//infile,form='unformatted',access='stream',status='old')

write(*,*) ' Enter the x and y grid dimensions:'
read(*,*) ntx,nty

call grabdata(ntx,nty)

!============================================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine grabdata(ntx,nty)
 !Reads data from unit 44 and processes

implicit none
integer:: ntx,nty,iopt,k,i
real:: t,rqa(nty,ntx)

read(44) t,rqa

write(*,*) ' Choose:'
write(*,*) ' (1) Constant  latitude cross section, or'
write(*,*) ' (2) Constant longitude cross section?'
read(*,*) iopt

write(*,*) ' Grid point index?'
read(*,*) k

open(55,file='cross.asc',status='replace')
if (iopt .eq. 1) then
  do i=1,ntx
    write(55,'(i5,1x,e14.7)') i,rqa(k,i)
  enddo
else
  do i=1,nty
    write(55,'(i5,1x,e14.7)') i,rqa(i,k)
  enddo
endif

end subroutine

!=========================================================

end program
