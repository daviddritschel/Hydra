program jet
! ----------------------------------------------------------------------
! |   This routine sets up the initial PV field in each layer for a    |
! |   baroclinic or an upper-layer jet.                                |
!-----------------------------------------------------------------------

 !Import constants and parameters:
use constants

implicit double precision (a-h,o-z)

 !Local arrays:
double precision:: bb(0:ny,0:nxm1),qbc(0:ny,0:nxm1)
double precision:: qq(0:ny,0:nxm1,nz)

 !-----------------------------------------------
write(*,*) ' Eastward or westward jet (1 or -1)?'
read(*,*) idir
 !Choose dq = +/-4*pi:
dq=sign(one,dble(idir))*four*pi
hdq=dq/two

write(*,*) ' Type of jet, (1) baroclinic or (2) upper-layer only?'
read(*,*) ityp

write(*,*) ' Width of the jet?'
read(*,*) w
hw=w/two

 !PV gradient in jet:
pvg=dq/w

write(*,*) ' We perturb the centreline of the jet by eps*sin(2*pi*x/L_x).'
write(*,*) ' Enter eps:'
read(*,*) eps

 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
  open(12,file='topo.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) dum,bb
  close(12)
else
  bb=zero
endif

if (ityp .eq. 1) then
   !Set up the baroclinic (mode 2) PV distribution in the array qbc:
  do ix=0,nxm1
    x=xmin+glx*dble(ix)
    dy=eps*sin(twopi*x/ellx)
    y1=-hw+dy
    y2= hw+dy
    do iy=0,ny
      y=ymin+gly*dble(iy)
      if (y .lt. y1) then
        qbc(iy,ix)=-hdq
      else if (y .gt. y2) then
        qbc(iy,ix)= hdq
      else
        qbc(iy,ix)=pvg*(y-dy)
      endif
    enddo
  enddo

   !Convert to PV in each layer, then add on beta*y and topography:
  do ix=0,nxm1
    do iy=0,ny
      qbt=beta*(ymin+gly*dble(iy))
      qq(iy,ix,1)=qbt+vect12*qbc(iy,ix)+bb(iy,ix)
      qq(iy,ix,2)=qbt+vect22*qbc(iy,ix)
    enddo
  enddo

else
   !Upper layer jet only; confine PV variations there:
  do ix=0,nxm1
    x=xmin+glx*dble(ix)
    dy=eps*sin(twopi*x/ellx)
    y1=-hw+dy
    y2= hw+dy
    do iy=0,ny
      y=ymin+gly*dble(iy)
      if (y .lt. y1) then
        qq(iy,ix,2)=-hdq
      else if (y .gt. y2) then
        qq(iy,ix,2)= hdq
      else
        qq(iy,ix,2)=pvg*(y-dy)
      endif
    enddo
  enddo

   !Add on beta*y and topography:
  do ix=0,nxm1
    do iy=0,ny
      qbt=beta*(ymin+gly*dble(iy))
      qq(iy,ix,1)=qbt+bb(iy,ix)
      qq(iy,ix,2)=qbt+qq(iy,ix,2)
    enddo
  enddo
endif

 !Write initial PV distribution to file:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
do iz=1,nz
  write(11,rec=iz) t,qq(:,:,iz)
enddo
close(11)

 !Write input data to file:
open(12,file='input_for_jet',status='unknown')
if (idir .gt. 0) then
  write(12,*) '  Eastward jet with dq = 4*pi'
  write(12,'(4x,f6.4,a)')  h2*pvg/beta,' ! h_2*dq/beta*w'
else
  write(12,*) '  Westward jet with dq = -4*pi'
  write(12,'(4x,f6.4,a)') -h1*pvg/beta,' ! -h_1*dq/beta*w'
endif
if (ityp .gt. 1) then
  write(12,*) '  Baroclinic jet.'
else
  write(12,*) '  Upper layer jet.'
endif
write(12,'(4x,f9.5,a)')   h1,' ! lower layer fractional thickness'
write(12,'(4x,f9.5,a)')   h2,' ! upper layer fractional thickness'
write(12,'(4x,f9.5,a)') beta,' ! planetary vorticity gradient, beta'
write(12,'(4x,f9.5,a)')    w,' ! width of the jet'
write(12,'(4x,f9.5,a)')  eps,' ! amplitude of centreline displacement'
close(12)

end program jet
