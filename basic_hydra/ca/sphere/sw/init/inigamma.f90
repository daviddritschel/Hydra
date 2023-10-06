program inigamma

!=====================================================
! Initialises the acceleration divergence given the
! PV, the height and the divergence fields.  
!=====================================================

 !Import spectral module:
use spectral

 !Declarations:
implicit none

double precision:: zz(ng,nt),dd(ng,nt),hh(ng,nt)
double precision:: zs(ng,nt),ds(ng,nt),hs(ng,nt)
double precision:: ud(ng,nt),uu(ng,nt),gg(ng,nt)
double precision:: wka(ng,nt),wkb(ng,nt)
double precision:: t
integer:: i,m

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

!----------------------------------------------------------------------
 !Read in the dimensionless height anomaly, hh:
open(11,file='hh_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,hh
close(11)

 !Read in the velocity divergence, dd:
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,dd
close(11)

 !Read in gridded PV, qs = (zeta+f)/(1+h), where zeta is the relative
 !vorticity and h is the dimensionless height anomaly:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)

 !Define relative vorticity:
do i=1,nt
  zz(:,i)=(one+hh(:,i))*zz(:,i)-cof
enddo

 !Create a spectral versions of hh, zz and dd:
zs=zz
call forfft(ng,nt,zs,trig,factors) 
ds=dd
call forfft(ng,nt,ds,trig,factors) 
hs=hh
call forfft(ng,nt,hs,trig,factors) 

 !Compute Laplacian of h:
call laplace(hs,gg) ! gg is returned in semi-spectral space
 !Bring gg back to physical space:
call revfft(ng,nt,gg,trig,factors) 

 !Invert Laplace's operator on the divergence (ds):
call laplinv(ds,wka,wkb)
 !Here the divergence potential is wka while d(wka)/dlat = wkb,
 !the divergent meridional velocity (all in semi-spectral space)

 !Compute divergent zonal velocity and store in ud:
call deriv(ng,nt,rk,wka,ud) 
do m=1,nt
  ud(:,m)=clati*ud(:,m)
enddo
 !ud = (1/r)*d(wka)/dlon where r = cos(lat)

 !Invert Laplace's operator on the relative vorticity (zs):
call laplinv(zs,wka,wkb)
 !Here the streamfunction is wka while d(wka)/dlat = wkb.

 !Complete calculation of zonal velocity, uu:
uu=ud-wkb
 !Bring uu back to physical space:
call revfft(ng,nt,uu,trig,factors) 

 !Add everything up to define acceleration divergence, gg:
do i=1,nt
  gg(:,i)=cof*zz(:,i)-bet*uu(:,i)-csq*gg(:,i)
enddo

 !Write data:
open(20,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) 0.d0,gg
close(20)

write(*,*) ' Acceleration divergence initialised in gg_init.r8'
end program inigamma
