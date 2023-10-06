program poten
!  -------------------------------------------------------------------------
!  |   Computes the magnetic tension from the inverse Laplacian of the     |
!  |   curl of the lorentz force, and the magnetic pressure from the       |
!  |   inverse Laplacian of minus the divergence of the lorentz force,     |
!  |   using data in jj.r4.                                                |
!  |                                                                       |
!  |   Output is to the unformatted, direct-access files "tt.r4" and       |
!  |   "pp.r4", which can be viewed with dataview just like jj.r4.         |
!  -------------------------------------------------------------------------

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

double precision:: pp(ny,nx),jj(ny,nx)
double precision:: jjx(ny,nx),jjy(ny,nx),aax(ny,nx),aay(ny,nx)
double precision:: ss(nx,ny),vtmp(nx,ny)
real:: qqr4(ny,nx),t
integer:: loop,iread,ix,iy,kx,ky

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open file containing current density field:
open(31,file='evolution/jj.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output files:
open(41,file='evolution/pp.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(42,file='evolution/tt.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) t,qqr4
  if (iread .ne. 0) exit 

   !Compute j & its derivatives:
  do ix=1,nx
    do iy=1,ny
      jj(iy,ix)=dble(qqr4(iy,ix))
    enddo
  enddo

  call ptospc(nx,ny,jj,ss,xfactors,yfactors,xtrig,ytrig)
   !Apply de-aliasing filter:
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=filt(kx,ky)*ss(kx,ky)
    enddo
  enddo

  call xderiv(nx,ny,hrkx,ss,vtmp)
  call spctop(nx,ny,vtmp,jjx,xfactors,yfactors,xtrig,ytrig)

  call yderiv(nx,ny,hrky,ss,vtmp)
  call spctop(nx,ny,vtmp,jjy,xfactors,yfactors,xtrig,ytrig)

   !Compute A from -grad^{-2}(j):
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=-green(kx,ky)*ss(kx,ky)
    enddo
  enddo

   !Get derivatives of A:
  call xderiv(nx,ny,hrkx,ss,vtmp)
  call spctop(nx,ny,vtmp,aax,xfactors,yfactors,xtrig,ytrig)

  call yderiv(nx,ny,hrky,ss,vtmp)
  call spctop(nx,ny,vtmp,aay,xfactors,yfactors,xtrig,ytrig)

   !Compute magnetic pressure (pp) and magnetic tension (re-use jj):
  do ix=1,nx
    do iy=1,ny
      pp(iy,ix)=jj(iy,ix)**2-aax(iy,ix)*jjx(iy,ix)-(b0+aay(iy,ix))*jjy(iy,ix)
      jj(iy,ix)=(b0+aay(iy,ix))*jjx(iy,ix)-aax(iy,ix)*jjy(iy,ix)
    enddo
  enddo

   !De-aliase and invert Laplacian to get magnetic pressure:
  call ptospc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=filt(kx,ky)*green(kx,ky)*ss(kx,ky)
    enddo
  enddo
  call spctop(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

   !Convert to real*4:
  qqr4=real(pp)
   !Write data:
  write(41,rec=loop) t,qqr4

   !De-aliase and invert Laplacian to get magnetic tension:
  call ptospc(nx,ny,jj,ss,xfactors,yfactors,xtrig,ytrig)
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=filt(kx,ky)*green(kx,ky)*ss(kx,ky)
    enddo
  enddo
  call spctop(nx,ny,ss,jj,xfactors,yfactors,xtrig,ytrig)

  write(*,'(a,f12.6)') ' t = ',t

   !Convert to real*4:
  qqr4=real(jj)
   !Write data:
  write(42,rec=loop) t,qqr4

enddo

 !Close files:
close(31)
close(41)
close(42)

write(*,*)
write(*,*) ' The magnetic pressure is ready in pp.r4 while the magnetic'
write(*,*) ' tension is ready in tt.r4'

return
end subroutine

 !End main program
end program
!=======================================================================
