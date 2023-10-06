program source
!  -------------------------------------------------------------------------
!  |   Computes the curl of the lorentz force from data in jj.r4           |
!  |                                                                       |
!  |   Output is to the unformatted, direct-access file "ss.r4" which      |
!  |   can be viewed with dataview just like qq.r4 and jj.r4               |
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

double precision:: jj(ny,nx),jjx(ny,nx),jjy(ny,nx),aax(ny,nx),aay(ny,nx)
double precision:: ss(nx,ny),vtmp(nx,ny),cmult,ssl1,ssmin,ssmax
real:: qqr4(ny,nx),t
integer:: loop,iread,ix,iy,kx,ky

write(*,*) ' The field values are capped at a multiple C of the L1 norm.'
write(*,*) ' Enter C:'
read(*,*) cmult

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open file containing current density field:
open(31,file='evolution/jj.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output file:
open(41,file='evolution/ss.r4',form='unformatted', &
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

   !Compute source term (re-use jj):
  do ix=1,nx
    do iy=1,ny
      jj(iy,ix)=(b0+aay(iy,ix))*jjx(iy,ix)-aax(iy,ix)*jjy(iy,ix)
    enddo
  enddo

   !Apply de-aliasing filter:
  call ptospc(nx,ny,jj,ss,xfactors,yfactors,xtrig,ytrig)
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=filt(kx,ky)*ss(kx,ky)
    enddo
  enddo
  call spctop(nx,ny,ss,jj,xfactors,yfactors,xtrig,ytrig)

   !Compute L1 norm and cap values:
  ssl1=zero
  do ix=1,nx
    do iy=1,ny
      ssl1=ssl1+abs(jj(iy,ix))
    enddo
  enddo
  ssl1=ssl1*dsumi

  write(*,'(a,f12.6,a,1p,e14.7)') ' t = ',t,'   ||S||_1 = ',ssl1

  ssmax=cmult*ssl1
  ssmin=-ssmax
  do ix=1,nx
    do iy=1,ny
      jj(iy,ix)=min(ssmax,max(ssmin,jj(iy,ix)))
    enddo
  enddo

   !Convert to real*4:
  qqr4=real(jj)

   !Write data:
  write(41,rec=loop) t,qqr4

enddo

 !Close files:
close(31)
close(41)

write(*,*)
write(*,*) ' The results are ready in ss.r4'

return
end subroutine

 !End main program
end program
!=======================================================================
