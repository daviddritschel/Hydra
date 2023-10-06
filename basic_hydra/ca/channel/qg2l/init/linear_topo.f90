program linear_topo
! ---------------------------------------------------------------------|
! |   This routine sets up an initial flow (U,0) which is uniform in   |
! |   each layer (as prescribed in the flow-setup script), apart from  |
! |   a small Rossby wave perturbation in the lower layer.  Also sets  |
! |   up linear topography which is compensated by the lower layer PV. |
!----------------------------------------------------------------------|

 !Import constants & parameters:
use constants

implicit double precision (a-h,o-z)

double precision:: qs(0:ny,0:nxm1)

if (.not. topogr) then
   write(*,*) ' This routine requires the topographic switch to be true!'
   write(*,*) ' *** exiting ***'
   stop
endif

write(*,*) ' We consider topography of the form H_b = sigma*y where sigma'
write(*,*) ' is the topographic slope.  The QG equations use the topographic'
write(*,*) ' PV, defined by q_b = f_0*H_b/H_1 where H_1 is the depth of the'
write(*,*) ' lower layer.  Hence, we specify the topographic "beta", defined'
write(*,*) ' by beta_b = f_0*sigma/H_1.'
write(*,*)
write(*,*) ' Enter beta_b/beta:'
read(*,*) betabnd

betab=beta*betabnd

 !Ask for a perturbation in the form of a pure Rossby wave:
write(*,*) ' We add a pure Rossby wave to the PV field in layer 1 of the form'
write(*,*) '        B*sin(k_n*x)*sin(l_m*y) '
write(*,*) ' where k_n = 2*pi*n/L_x and l_m=pi*m/L_y.'
write(*,*) ' Enter B, n & m:'
read(*,*) bross,nross,mross

akx=twopi*dble(nross)/ellx
aky=   pi*dble(mross)/elly

 !---------------------------------------------------------------------------
 !Write initial PV distribution to a file 
 !(NOTE: the centre of the channel is y = 0):
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
 !Layer 1:
pvg1=epsilon*beta-betab
do ix=0,nxm1
  xg=xmin+glx*dble(ix)
  do iy=0,ny
    yg=ymin+gly*dble(iy)
    qs(iy,ix)=pvg1*yg + bross*sin(akx*xg)*sin(aky*yg)
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
 !Write topography to a file (NOTE: this must have zero domain mean):
do ix=0,nxm1
  xg=xmin+glx*dble(ix)
  do iy=0,ny
    yg=ymin+gly*dble(iy)
    qs(iy,ix)=betab*yg
  enddo
enddo

open(11,file='topo.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qs
close(11)
   
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

end program linear_topo
