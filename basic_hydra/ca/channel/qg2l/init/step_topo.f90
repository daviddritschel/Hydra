program step_topo
! --------------------------------------------------------------------
! |   This routine sets up an initial flow at rest; it is used to    |
! |   set up a step typography.                                      |
!---------------------------------------------------------------------

 !Import constants & parameters:
use constants

implicit double precision (a-h,o-z)

double precision:: qs(0:ny,0:nxm1)

 !---------------------------------------------------------------------------
 !Write initial PV distribution to a file 
 !(NOTE: the centre of the channel is y = 0):
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
 !Layer 1:
pvg1=epsilon*beta
do ix=0,nxm1
  xg=xmin+glx*dble(ix)
  do iy=0,ny
    yg=ymin+gly*dble(iy)
    qs(iy,ix)=pvg1*yg
  enddo
enddo
write(11,rec=1) zero,qs
 !Layer 2:
pvg2=lambda*beta
do ix=0,nxm1
  do iy=0,ny
    qs(iy,ix)=pvg2*(ymin+gly*dble(iy))
  enddo
enddo
write(11,rec=2) zero,qs
close(11)

 !---------------------------------------------------------------------------
 !Deal with step topography if present:
if (topogr) then
   !Add topography (NOTE: this must have zero domain mean):
  write(*,*) ' We consider step topography H_b of the form '
  write(*,*) '   f_0*H_b/H_1 = A for yc-D < y < yc+D where '
  write(*,*) ' H_1 is the depth of the lowest layer, '
  write(*,*) ' y_c = (y_min+y_max)/2 and D the '
  write(*,*) ' latitudinal width of the step. '
  write(*,*)
  write(*,*) ' Enter A & D:'
  read(*,*) topamp,displa

  write(*,*) ' Enter amplitude of random noise:'
  read(*,*) dtopamp

  hdispla=displa/two

  qsbar=zero
  do iy=0,ny
    yd=ymin+gly*dble(iy)-ycen
    if (abs(yd) .lt. hdispla) then
      do ix=0,nxm1
        qs(iy,ix)=topamp
      enddo
      qsbar=qsbar+topamp
    else
      do ix=0,nxm1
        qs(iy,ix)=zero
      enddo
    endif
  enddo
  qsbar=qsbar/dble(ny)

   !Remove mean and add random noise:
  do ix=0,nxm1
    do iy=0,ny
      qs(iy,ix)=qs(iy,ix)-qsbar+dtopamp*(two*rand(0)-one)
    enddo
  enddo

   !Write topography to a file:
  open(11,file='topo.r8',form='unformatted', &
      & access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,qs
  close(11)
endif
   
 !---------------------------------------------------------------------------
if (heating) then
   !Relax back to slanting interfaces:

   !Obtain mean velocity in each layer given fraction PV gradient epsilon

   !Middle interface (called "1"):
  slope=-h1h2kdbarsq*(u1lmean-alpha*u2lmean)
   !Note: dd1(iy,ix)=h1h2kdbarsq*(pp(iy,ix,1)-alpha*pp(iy,ix,2)) in evolution.f90
  do ix=0,nxm1
    do iy=0,ny
      qs(iy,ix)=slope*(ymin+gly*dble(iy))
    enddo
  enddo

   !Write displacement to a file:
  open(11,file='disp1eq.r8',form='unformatted', &
      & access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,qs
  close(11)

  if (.not. barot) then
   !Upper interface (called "2"):
    slope=-h1h2ackdbarsq*u2lmean
   !Note: dd2(iy,ix)=h1h2ackdbarsq*pp(iy,ix,2) in evolution.f90
    do ix=0,nxm1
      do iy=0,ny
        qs(iy,ix)=slope*(ymin+gly*dble(iy))
      enddo
    enddo

   !Write displacement to a file:
    open(12,file='disp2eq.r8',form='unformatted', &
        & access='direct',status='replace',recl=2*nbytes)
    write(12,rec=1) zero,qs
    close(12)
  endif
endif

end program
