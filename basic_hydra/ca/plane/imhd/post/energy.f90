program energy
!  ----------------------------------------------------------------------
!  |   Computes |u^2| and |B^2| from data in qq.r4 and aa.r4            |
!  |                                                                    |
!  |   Output is to the unformatted direct-access files "usq.r4" and    |
!  |   "bsq.r4", which can be viewed with dv just like qq.r4, etc.      |
!  ----------------------------------------------------------------------

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use spectral

implicit none

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

double precision:: uu(ny,nx),vv(ny,nx),qq(ny,nx)
double precision:: aa(ny,nx),aax(ny,nx),aay(ny,nx)
double precision:: usq(ny,nx),bsq(ny,nx)
double precision:: wka(nx,ny),wkb(nx,ny)
real:: qqr4(ny,nx),t
integer:: loop,iread

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open file containing q & A:
open(31,file='evolution/qq.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(32,file='evolution/aa.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output files:
open(41,file='evolution/usq.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(42,file='evolution/bsq.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) t,qqr4
  if (iread .ne. 0) exit

  qq=dble(qqr4)
  call ptospc(nx,ny,qq,wka,xfactors,yfactors,xtrig,ytrig)
  call main_invert(wka,uu,vv,wkb)
  usq=uu**2+vv**2
  
  read(32,rec=loop) t,qqr4
  aa=dble(qqr4)

  call ptospc(nx,ny,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call gradient(wka,aax,aay)
  bsq=aax**2+aay**2

   !Write data:
  write(41,rec=loop) t,real(usq)
  write(42,rec=loop) t,real(bsq)
enddo

 !Close files:
close(31)
close(32)
close(41)
close(42)

write(*,*)
write(*,*) ' |u|^2 & |B^2| are available in evolution/usq.r4 & bsq.r4.'
write(*,*)

return
end subroutine

 !End main program
end program energy
!=======================================================================
