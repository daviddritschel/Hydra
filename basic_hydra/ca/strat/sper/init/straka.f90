program straka

use constants

! This routine sets up initial buoyancy and vorticity fields for the
! Straka test case (see EPIC-2D paper).

implicit none

double precision,parameter:: T0=300.d0 !Reference temperature (degrees K)
double precision,parameter:: g=9.80665 !Acceleration due to gravity (m/s^2)
double precision:: zz(0:ny,0:nxm1),bb(0:ny,0:nxm1)
double precision:: x0,y0,a0,b0,dTmax,bfac
double precision:: x,y,r
integer:: ix,iy

!---------------------------------------------------------------------
write(*,'(a,2(f6.2),a)') &
     ' We consider a domain of height and length: ', &
     elly/1000.d0,ellx/1000.d0,'km.'
write(*,'(a,2(i6))') ' The y and x grid resolution is: ',ny,nx
write(*,*)

write(*,*) ' Enter the centre (x_o,y_o) of the elliptical cold bubble (in km):'
read(*,*) x0,y0
x0=1000.d0*x0
y0=1000.d0*y0

write(*,*) ' Enter the semi-axes (a_o,b_o) of the ellipse (in km):'
read(*,*) a0,b0
a0=1000.d0*a0
b0=1000.d0*b0

write(*,*) ' Enter the temperature perturbation Delta T_max (in degrees K):'
read(*,*) dTmax

bfac=-g*dTmax/(two*T0)

 !Set up buoyancy distribution:
do ix=0,nxm1
  x=(xmin+glx*dble(ix)-x0)/a0
  do iy=0,ny
    y=(ymin+gly*dble(iy)-y0)/b0
    r=sqrt(x**2+y**2)
    if (r < one) then
      bb(iy,ix)=bfac*(cos(pi*r)+one)
    else
      bb(iy,ix)=zero
    endif
    zz(iy,ix)=zero
  enddo
enddo

 !Write vorticity distribution to file:
open(20,file='zz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

 !Write buoyancy distribution to file:
open(20,file='bb_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,bb
close(20)

 !Write input data to file:
x0=x0/1000.d0
y0=y0/1000.d0
a0=a0/1000.d0
b0=b0/1000.d0
open(12,file='input_for_straka',status='unknown')
write(12,'(4x,f5.1,a)') x0,' ! centre in x of the elliptical cold bubble (km)'
write(12,'(4x,f5.1,a)') y0,' ! centre in y "   "      "       "     "     "  '
write(12,'(4x,f5.1,a)') a0,' ! semi-major axis of the cold bubble (km)'
write(12,'(4x,f5.1,a)') b0,' ! semi-minor axis "   "   "     "     "  '
write(12,'(4x,f5.1,a)') dTmax,' ! maximum temperature perturbation'
close(12)
      
end program
