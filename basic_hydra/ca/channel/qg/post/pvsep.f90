program pvsep
!  -------------------------------------------------------------------------
!  |   Divides the PV anomaly, q - beta*y, into a symmetric and an         |
!  |   antisymmetric part, writes these for imaging as qqa & qqs.r4        |
!  |   respectively, and computes the associated enstrophy, writing        |
!  |   this (and the total enstrophy) to enstrophy.asc.                    |
!  -------------------------------------------------------------------------

 !Import constants and parameters:
use constants

implicit none

integer,parameter:: nyh=ny/2,nhbytes=4*(nx*(nyh+1)+1)

!---------------------------------------------------------
 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(0:ny,0:nxm1),qa(0:nyh,0:nxm1),qs(0:nyh,0:nxm1)

 !Diagnostic quantities:
double precision:: bety(0:ny),enst,ensa,enss
real:: qqr4(0:ny,0:nxm1)
real:: tr4

integer:: loop,iread,ix,iy

!---------------------------------------------------------------
 !Open file containing PV field:
open(31,file='evolution/qq.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output files:
open(22,file='evolution/enstrophy.asc',status='replace')

open(41,file='evolution/qa.r4',form='unformatted', &
        access='direct',status='replace',recl=nhbytes)

open(42,file='evolution/qs.r4',form='unformatted', &
        access='direct',status='replace',recl=nhbytes)

 !Define beta*y:
do iy=0,ny
  bety(iy)=beta*(ymin+gly*dble(iy))
enddo

!---------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit 
  write(*,'(a,f12.5)') ' Processing t = ',tr4

   !Convert PV to double precision as qq:
  qq=dble(qqr4)

   !Remove beta*y:
  do ix=0,nxm1
    do iy=0,ny
      qq(iy,ix)=qq(iy,ix)-bety(iy)
    enddo
  enddo

   !Define qa & qs:
  do ix=0,nxm1
    do iy=0,nyh
      qa(iy,ix)=f12*(qq(nyh+iy,ix)-qq(nyh-iy,ix))
      qs(iy,ix)=f12*(qq(nyh+iy,ix)+qq(nyh-iy,ix))
    enddo
  enddo

   !Calculate various enstrophies:
  enst=f12*garea*sum(qq**2)
  ensa=garea*sum(qa**2)
  enss=garea*sum(qs**2)

   !Write data:
  write(22,'(1x,f13.5,3(1x,e14.7))') tr4,enst,ensa,enss

  write(41,rec=loop) tr4,real(qa)
  write(42,rec=loop) tr4,real(qs)

enddo

 !Close output files:
close(22)
close(41)
close(42)

write(*,*)
write(*,*) ' The results are ready in evolution/qa.r4, evolution/qs.r4,'
write(*,*) ' and evolution/enstrophy.asc'

return
end subroutine diagnose

 !End main program
end program pvsep
!=======================================================================
