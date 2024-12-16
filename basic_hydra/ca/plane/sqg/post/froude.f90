program froude
!  -------------------------------------------------------------------------
!  |      Computes Fr = Ri^{-1/2} = max{|grad{q}|} from data in bb.r4      |
!  |      Creates froude.asc                                               |
!  -------------------------------------------------------------------------

 !Import constants, parameters and arrays needed for FFTs:
use spectral

implicit none

!---------------------------------------------------------
 !Initialise FFTs:
call init_spectral

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(ny,nx),qx(ny,nx),qy(ny,nx)

 !Spectral field:
double precision:: wka(nx,ny)

 !Work quantities:
real:: qqr4(ny,nx)
real:: t
double precision:: frmax

integer:: loop,iread

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open file containing the scaled buoyancy field:
open(31,file='bb.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output file:
open(41,file='froude.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=1
do
   iread=0
   read(31,rec=loop,iostat=iread) t,qqr4
   if (iread .ne. 0) exit 
   write(*,'(a,f12.5)') ' Processing t = ',t

   qq=dble(qqr4)

    !Convert to spectral space:
   call ptospc(nx,ny,qq,wka,xfactors,yfactors,xtrig,ytrig)

    !Compute gradient:
   call gradient(wka,qx,qy)

    !Compute Froude number:
   qq=sqrt(qx**2+qy**2)
   frmax=maxval(qq)

    !Write diagnostic data:
   write(41,'(1x,f12.5,1x,1p,e14.7)') t,frmax

   loop=loop+1
enddo

 !Close files:
close(31)
close(41)

write(*,*)
write(*,*) ' The diagnostics are in froude.asc'

return
end subroutine diagnose

 !End main program
end program froude
!=======================================================================
