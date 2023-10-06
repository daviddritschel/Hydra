program staircase

! Sets up a multiple Bickley jets, superposed with a small non-zonal
! perturbation.
! IC adapted from one suggested by Inna Polichtchouk & James Cho in 
! February 2013

use constants 

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: zz(ng,nt)
double precision:: obar(ng),zbar(ng),clat(ng),slat(ng),cof(ng)

!-----------------------------------------------------------------------
write(*,*) ' We consider a planet of radius 1 rotating with a period of'
write(*,*) ' one "day", i.e. at a rotation rate of 2*pi.'
write(*,*)
write(*,*) ' We start with a superposition of N Bickley jets equally spaced'
write(*,*) ' in latitude, with maximum angular velocity 2*pi*eps.'
write(*,*) ' Enter N & eps:'
read(*,*) nj,eps
write(*,*)
write(*,*) ' We add a vorticity perturbation of the form 2*pi*A*exp(-4*d^2/(4-d^2)),'
write(*,*) ' where d is the chord distance from a point at latitude,longitude (p,0).'
write(*,*) ' *** Note, the global vorticity perturbation mean is removed.'
write(*,*) ' Enter A & p (degrees):'
read(*,*) a,pz
pz=pz*pi/180.d0

!------------------------------------------------------------
hdl=f12*dl
 !dl: the latitude & longitude grid spacing
!------------------------------------------------------------
 !Define cos, sin and Coriolis frequency:
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
  cof(j)=fpole*slat(j)
enddo

 !For computing global means:
rsum=f1112*(clat(1)+clat(ng))
do j=2,ng-1
  rsum=rsum+clat(j)
enddo
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.

!----------------------------------------------------------
 !Define zonal angular velocity and absolute vorticity
 !from the superposition of jets:
do j=1,ng
  zbar(j)=cof(j)
enddo

omega0=twopi*eps
dphic=pi/dble(nj)
b=dphic*f14
binv=one/b

do k=1,nj
  phi0=dphic*(dble(k)-f12)-hpi
  do j=1,ng
    phi=dl*(dble(j)-f12)-hpi
    chlat=cosh((phi-phi0)*binv)
    shlat=sinh((phi-phi0)*binv)
    dobar=omega0/chlat**2
    obar(j)=obar(j)+dobar
    zbar(j)=zbar(j)+two*dobar*(slat(j)+binv*clat(j)*shlat/chlat)
  enddo
enddo

 !Remove global mean angular velocity:
osum=f1112*(clat(1)*obar(1)+clat(ng)*obar(ng))
do j=2,ng-1
  osum=osum+clat(j)*obar(j)
enddo
oavg=osum/rsum
do j=1,ng
  obar(j)=obar(j)-oavg
  zbar(j)=zbar(j)-two*oavg*slat(j)
enddo

!-------------------------------------------------------------------------
 !Record zonal profiles vs z = sin(latitude) for diagnostic purposes:
open(22,file='uz_zonal.asc',status='unknown')
do j=1,ng
  write(22,'(3(1x,f15.9))') slat(j),clat(j)*obar(j),zbar(j)-cof(j)
enddo
close(22)

!-------------------------------------------------------------------------
 !Construct vorticity perturbation:
zamp=twopi*a
cpz=cos(pz)
spz=sin(pz)
do i=1,nt
  rlon=dl*dble(i-1)-pi
  clon=cos(rlon)
  slon=sin(rlon)
  do j=1,ng
    zp=cpz*clat(j)*clon+spz*slat(j)
    arg=four*(one-zp)/(abs(one+zp)+small)
    if (arg .lt. 33.d0) then
      zz(j,i)=zamp*exp(-arg)
    else
       !exp(-arg) < machine double precision:
      zz(j,i)=zero
    endif
  enddo
enddo

 !Remove global mean:
zsum=zero
do i=1,nt
  zsum=zsum+f1112*(zz(1,i)*clat(1)+zz(ng,i)*clat(ng))
  do j=2,ngm1
    zsum=zsum+zz(j,i)*clat(j)
  enddo
enddo
zsum=zsum*dsumi

do i=1,nt
  do j=1,ng
    zz(j,i)=zz(j,i)-zsum
  enddo
enddo

!-------------------------------------------------------------------------
 !Define absolute vorticity:
do i=1,nt
  do j=1,ng
    zz(j,i)=zbar(j)+zz(j,i)
  enddo
enddo

 !Write initial absolute vorticity field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

end program
