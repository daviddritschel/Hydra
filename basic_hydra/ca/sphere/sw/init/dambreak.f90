program dambreak
!============================================================
! Uses a height anomaly field proportional to a tanh function
! centred on a given latitude.
!============================================================

use constants 

 !Declarations:
implicit none

double precision:: hh(ng,nt),qq(ng,nt)
double precision:: slat(ng),clat(ng)
double precision:: eps,amp,phi0,b,binv,phi,fm,acmlon
double precision:: rsum,rsumi,dsumi,vsum
integer:: i,j,m

!------------------------------------------------------------------------------
write(*,*) ' We consider an initial height anomaly field of the form'
write(*,*) '          h = eps*tanh((phi-phi_0)/b) + C'
write(*,*) ' where C is chosen so that h has zero global mean.'

write(*,*) ' Enter eps:'
read(*,*) eps

write(*,*) ' Enter phi_0 in degrees: '
read(*,*) phi0
phi0=phi0*pi/180.d0

write(*,*) ' Enter b in degrees: '
read(*,*) b
b=b*pi/180.d0

write(*,*)
write(*,*) ' To make the flow azonal, we displace the latitudes by'
write(*,*) ' A*cos(phi)*cos(m*lambda).'
write(*,*)
write(*,*) ' Enter the amplitude A (degrees): '
read(*,*) amp
amp=amp*pi/180.d0

write(*,*) ' Enter the wavenumber m (integer): '
read(*,*) m

!------------------------------------------------------------------------------
 !Define cos and sin(latitude):
do j=1,ng
  phi=dl*(dble(j)-f12)-hpi
  slat(j)=sin(phi)
  clat(j)=cos(phi)
enddo

 !For removing global averages:
rsum=f1112*(clat(1)+clat(ng))
do j=2,ngm1
  rsum=rsum+clat(j)
enddo
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.
rsumi=one/rsum
dsumi=rsumi/dble(nt)

 !Write initial divergence field:
hh=zero
open(20,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

!------------------------------------------------------------------------------
 !Compute initial height field:
fm=dble(m)
binv=one/b
do i=1,nt
  acmlon=amp*cos(fm*(dl*dble(i-1)-pi))
  do j=1,ng
    phi=dl*(dble(j)-f12)-hpi-clat(j)*acmlon
    hh(j,i)=eps*tanh(binv*(phi-phi0))
  enddo
enddo

 !Remove global mean height anomaly:
vsum=zero
do i=1,nt
  vsum=vsum+f1112*(clat(1)*hh(1,i)+clat(ng)*hh(ng,i))
  do j=2,ngm1
    vsum=vsum+clat(j)*hh(j,i)
  enddo
enddo
vsum=vsum*dsumi
do i=1,nt
  do j=1,ng
    hh(j,i)=hh(j,i)-vsum
  enddo
enddo

 !Write initial height field:
open(20,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

!Compute PV field for an initially resting state:
do i=1,nt
  do j=1,ng
    qq(j,i)=fpole*slat(j)/(one+hh(j,i))
  enddo
enddo

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program dambreak
