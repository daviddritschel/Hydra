program rossby
! Test case 6 Rossby-Haurwitz Wave [Williamson et al JCP 102 (1992)]

implicit none
integer:: ng

! Specify resolution:
write(*,*) ' Number of latitudinal divisions?'
read(*,*) ng

call gendata(ng)


contains 

!=======================================================================

subroutine gendata(ng)

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0
double precision,parameter:: hpi=pi/2.d0,twopi=2.d0*pi,pif=pi/180.d0
double precision,parameter:: g=9.80616d0       ! gravity (m/s^2)
double precision,parameter:: omega=7.292d-5    ! Omega_e (1/s)
double precision,parameter:: tday=twopi/omega  ! Length of day (s)
double precision,parameter:: ae=6371220.d0     ! Earth radius (m)

! f = dimensionless Coriolis frequency below
double precision,parameter:: ve=ae/tday, f=4.d0*pi

! Define constants used in the Rossby-Haurwitz test:
double precision,parameter:: w=7.848d-6*tday, aK=w, h0=8000.d0
integer,parameter:: mR=4

double precision:: clat(ng),slat(ng)
double precision:: hbar(ng),zbar(ng),dhdl(ng)
double precision:: A(ng),B(ng),C(ng),D(ng),E(ng)
double precision:: clon(2*ng),slon(2*ng),cRlon(2*ng),c2Rlon(2*ng)
double precision:: hh(ng,2*ng),dd(ng,2*ng),qq(ng,2*ng)

! Number of longitudinal divisions:
nt=2*ng

ngridp=ng*nt
nbytes=4*(ngridp+1)

!--------------------------------------------------------------
! Define perturbation used for the depth field:
write(*,*) ' We consider a perturbation of the form A*f_p*z_r'
write(*,*) ' to the vorticity field, where z_r is a rotated z'
write(*,*) ' coordinate with the NP at (lon0,lat0) and f_p is'
write(*,*) ' the polar value of the Coriolis frequency.'
write(*,*) ' Enter A (suggest 0.025):'
read(*,*) apert
write(*,*) ' Enter lon0 in degrees (suggest 50):'
read(*,*) rlonp
write(*,*) ' Enter lat0 in degrees (suggest 40):'
read(*,*) rlatp

dlon=twopi/dble(nt)
! dlon: the maximum grid length (at the equator)

dlat=pi/dble(ng)
! dlat: the grid spacing in latitude

!------------------------------------------------------------
! Constants:
aKs=aK**2
R=dble(mR)
R1=R+1.d0
R2=R+2.d0
Rs=R*R
wa=twopi+w
ca=w*(2.d0*twopi+w)/2.d0
caR=2.d0*Rs - R - 2.d0      
cb=2.d0*aK*wa/(R1*R2)
cbR=Rs+2.d0*R+2.d0
cc=0.25d0*aKs
ceR=Rs+3.d0*R+2.d0

!------------------------------------------------------------
! Define cos and sin(latitude):
do j=1,ng
  rlat=dlat*(dble(j)-0.5d0)-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
enddo

! Define cos(R*longitude) & cos(2*R*longitude):
do i=1,nt
  rlon=dlon*dble(i)-pi
  slon(i)=sin(rlon)
  clon(i)=cos(rlon)
  cRlon(i)=cos(R*rlon)
  c2Rlon(i)=cos(2.d0*R*rlon)
enddo

!----------------------------------------------------------
! Compute fields:
top=0.d0
bot=0.d0
do j=1,ng
  c2lat=clat(j)**2
  cRlat=clat(j)**R
  cfac=cc*c2lat**R
  A(j)=ca*c2lat + cfac*(R1*c2lat + caR - 2.d0*Rs/c2lat)
  top=top+A(j)*clat(j)
  bot=bot+clat(j)
  B(j)=cb*cRlat*(cbR - R1**2*c2lat)
  C(j)=cfac*(R1*c2lat - R2)
  D(j)=2.d0*wa*slat(j)
  E(j)=aK*slat(j)*cRlat*ceR
enddo

Abar=top/bot
havg = h0 + Abar*ve**2/g
write(*,'(a,f10.3)') ' Mean depth (m) = ',havg
csq=g*havg

! Compute Rossby Radius:
rrtmp=sqrt(csq)/f
rrad=rrtmp/ve

! Write basic parameters to a file:
open(10,file='basic.dat',status='unknown')
write(10,'(1x,f14.11,1x,f14.10)') f,rrad
close(10)

amp=ve**2/csq

r0=cos(pif*rlatp)
x0=r0*cos(pif*rlonp)
y0=r0*sin(pif*rlonp)
z0=sin(pif*rlatp)

do i=1,nt
  do j=1,ng
    x=clat(j)*clon(i)
    y=clat(j)*slon(i)
    z=slat(j)
    hh(j,i)=amp*(A(j)-Abar+B(j)*cRlon(i)+C(j)*c2Rlon(i))+apert*(x*x0+y*y0+z*z0)
    dd(j,i)=0.d0
    qq(j,i)=(D(j)-E(j)*cRlon(i))/(1.d0+hh(j,i))
  enddo
enddo

 !Write initial height field:
open(20,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) 0.d0,hh
close(20)

 !Write initial divergence field:
open(20,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) 0.d0,dd
close(20)

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) 0.d0,qq
close(20)

return 
end subroutine gendata

 !End main program
end program rossby
