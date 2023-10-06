program bickley

! Sets up a perturbed bickley jet just off the equator.
! IC devised by Inna Polichtchouk & James Cho in February 2013
! Revised source to f90 by S King @ St Andrews March 2013

use spectral

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: zz(ng,nt),zzp(ng,nt)
double precision:: zbar(ng),clat(ng),slat(ng)
double precision:: lambda

!--------------------------------------------------------------
write(*,*) ' We consider a planet of radius 1 rotating with'
write(*,*) ' a period of one "day".'
write(*,*)
write(*,*) ' We start with a Bickley jet'
write(*,*) '  u = u_0 sech^2((phi-phi_0)/b).'
write(*,*) ' Enter u_0/2*pi, phi_0 and b:'
read(*,*) u0nd,phi0,b
u0=u0nd*twopi

write(*,*)
write(*,*) ' We add a vorticity perturbation of the form'
write(*,*) '  z_0 exp(-(lambda/alpha)^2-((phi-phi_c)/beta)^2).'
write(*,*) ' Enter z_0/(2*Omega), phi_c, alpha and beta:'
read(*,*) z0nd,phic,alpha,beta
z0=fpole*z0nd

write(*,*) ' A further random vorticity perturbation is added.'
write(*,*) ' Enter the rms vorticity perturbation:'
read(*,*) eps
if (eps .gt. zero) then
   !Initialise spectral module:
  call init_spectral
   !Initialize random # generator:
  do i=1,iseed
    uni=rand(0)
  enddo
   !generate random relative vorticity field zzp with this rms value:
  call ranspec(zzp,eps)
   !Return zzp to physical space:
  call revfft(ng,nt,zzp,trig,factors)
endif

!------------------------------------------------------------
hdl=f12*dl
 !dl: the latitude & longitude grid spacing
!------------------------------------------------------------
 !Define cos, sin and tan(latitude):
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
enddo

!----------------------------------------------------------
 !Define zonal velocity and absolute vorticity:
q0=two*u0/b
do j=1,ng
  phi=dl*(dble(j)-f12)-hpi
  chlat=cosh((phi-phi0)/b)
  ubar=u0/chlat**2
  shlat=sinh((phi-phi0)/b)
  zbar(j)=q0*shlat/chlat**3+(ubar/clat(j)+fpole)*slat(j)
enddo

 !Define perturbation relative vorticity and remove global average:
do i=1,nt
  lambda=dl*dble(i-1)-pi
  do j=1,ng
    phi=dl*(dble(j)-f12)-hpi
    zz(j,i)=z0*exp(-(lambda/alpha)**2-((phi-phic)/beta)**2)
  enddo
enddo

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

 !Add on zonal absolute vorticity:
do i=1,nt
  do j=1,ng
    zz(j,i)=zbar(j)+zz(j,i)
  enddo
enddo

if (eps .gt. zero) then
   !Add on random vorticity perturbation:
  do i=1,nt
    do j=1,ng
      zz(j,i)=zz(j,i)+zzp(j,i)
    enddo
  enddo
endif

 !Write initial absolute vorticity field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

end program
