program cofm

! Program to calculate centre of mass of a billow
! Written on 20th Jun by sek and jm @ St Andrews

use common

implicit none

 !Main variable declarations:
integer:: loop,iread,imax,iy,ix
real:: tr4,bbr4(0:ny,0:nx)
double precision:: by(0:ny,0:nx),y(0:ny),bint,byint

 !Open input and output files
open(10,file='centres.asc',status='unknown')
open(12,file='bb.r4',form='unformatted',access='direct',status='old',recl=nbytes)

 !Define y:
do iy=0,ny
  y(iy)=gly*dble(iy)
enddo

loop=0
iread=0
 !Read each frame of field in turn and process:
do 
  loop=loop+1
  iread=0
  read(12,rec=loop,iostat=iread) tr4,bbr4
   !Upon any read error exit loop and end program: 
  if (iread .ne. 0) exit 
  do iy=0,ny
    do ix=0,nx
      by(iy,ix)=dble(bbr4(iy,ix))*y(iy)
      bb(iy,ix)=dble(bbr4(iy,ix))
    enddo
  enddo
  call doubleint(by,byint)
  call doubleint(bb,bint)
  byint=byint/(bint+small)
  write(10,'(2(f13.9))') tr4,byint
enddo

close(12)
close(10)


contains

!===========================================================================

subroutine doubleint(qq,val)

! Computes the double integral of a field qq and returns the result in val:
!          ql1 = int_xmin^xmax{int_ymin^ymax{qq dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(0:ny,0:nx),val
!---------------------------------------------------------------
 
!Use trapezoidal rule in both directions:
val=zero
do ix=1,nxm1
  val=val+qq(0,ix)+qq(ny,ix)
enddo
do iy=1,nym1
  val=val+qq(iy,0)+qq(iy,nx)
enddo
val=f12*val+f14*(qq(0,0)+qq(ny,0)+qq(0,nx)+qq(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    val=val+qq(iy,ix)
  enddo
enddo

val=garea*val
 !Note: garea is the grid box area, glx*gly
return 
end subroutine
!---------------------------------------------------------------------------------

end program
