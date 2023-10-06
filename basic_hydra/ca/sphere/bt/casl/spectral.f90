module spectral

use constants
use stafft
use deriv1d

implicit none

 !Common arrays, constants:

 !FFT commons:
integer:: factors(5)
double precision:: trig(2*nt),rk(nt)
double precision:: wave(nt)

 !Tridiagonal commons:
double precision:: d2lon(ng,nt)
double precision:: sigm(2),psig(2)
double precision:: amla,apla,pc12,pc24
integer:: ksig(nt)

double precision:: cmla(ng,nt),dmla(ng)
double precision:: ala0(ng,nt),bla0(ng,nt)
double precision:: cla0(ng,nt),dla0(ng,nt)
double precision:: ela1(ng,nt),ela2(ng,nt)
double precision:: fla1(ng,nt),fla2(ng,nt)
double precision:: etd0(ng),htd0(ng)
double precision:: rsumi,dsumi,dl24

 !Latitudinal arrays:
double precision:: rho(ng),rhoi(ng),bet(ng),cof(ng)
double precision:: tau(ng),tsqi(ng),tdr(ng),rdt(ng)

contains

!===========================
subroutine init_spectral

 !Declarations:
implicit none

 !Local array:
double precision:: tlat(ng),rtisq(ng)
double precision:: aspsqm1,rlat,slat,dlt,wm,fac,rsum
double precision:: sm,pm,a0,b0,c0,d0,deti,cpla,dpla
integer:: j,m,k

!------------------------------------------------------------------------
 !Latitude functions:
aspsqm1=asp**2-one
do j=1,ng
  rlat=(dble(j)-f12)*dl-hpi
   !sin(lat):
  slat=sin(rlat)
   !rho = cos(lat):
  rho(j)=cos(rlat)
   !1/rho:
  rhoi(j)=one/rho(j)
   !1/tau^2:
  tsqi(j)=one+aspsqm1*rho(j)**2
   !tau:
  tau(j)=one/sqrt(tsqi(j))
   !rho/tau:
  rdt(j)=rho(j)/tau(j)
   !tau/rho:
  tdr(j)=one/rdt(j)
   !tau'/tau:
  dlt=aspsqm1*tau(j)**2*slat*rho(j)
   !-(rho*tau)'/(rho*tau) = tau^2*tan(lat):
  tlat(j)=slat*tau(j)*tdr(j)
   !1/(rho*tau)^2:
  rtisq(j)=rhoi(j)**2+aspsqm1
   !Coriolis frequency, f = -2*Omega*rho'*tau = 2*Omega*sin(lat)*tau:
  cof(j)=fpole*slat*tau(j)
   !beta = df/dlat = 2*Omega*cos(lat)*tau^3*b^2:
  bet(j)=fpole*rho(j)*tau(j)+cof(j)*dlt
enddo

!----------------------
 !Initialise the FFT module:
call initfft(nt,factors,trig)
 !Initialise the spectral derivative module:
call init_deriv(nt,twopi,rk)

 !Define longitudinal wavenumbers:
wave(1)=zero
do m=2,ng
  wm=rk(2*(m-1))
  wave(     m)=wm
  wave(ntp2-m)=wm
enddo
wave(ng+1)=rk(nt)

 !m^2/(rho^2*tau^2) factor which appears in Laplace's operator:
do m=1,nt
  fac=wave(m)**2
  do j=1,ng
    d2lon(j,m)=fac*rtisq(j)
  enddo
enddo

!----------------------------------------------------------------
 !Initialise tridiagonal coefficients:

!-------------------------------------------------------------
 !Sign changes for longitudinal wavenumbers in inversion:

 !Odd physical wavenumbers (1, 3, ...):
do m=2,nt,2
  ksig(m)=1
enddo
sigm(1)=-one
psig(1)=0.9d0

 !Even physical wavenumbers (0, 2, ...):
do m=1,ntm1,2
  ksig(m)=2
enddo
sigm(2)= one
psig(2)=1.1d0

 !Constants used immediately below:
pc12=1.2d0/dl**2
pc24=two*pc12
dl24=dl/24.d0
amla=0.6d0/dl
apla=-amla

!-----------------------------------------------------
 !Initialise Laplace block-tridiagonal inversion:

 !zonally-symmetric part (physical wavenumber 0):
htd0(1)=one/f74
etd0(1)=-f14*htd0(1)

do j=2,ngm1
  htd0(j)=one/(f32+f14*etd0(j-1))
  etd0(j)=-f14*htd0(j)
enddo

htd0(ng)=one/(f74+f14*etd0(ngm1))

 !Used in integrals over mu where dmu = rho/tau:
rsum=f1112*(rdt(1)+rdt(ng))
do j=2,ngm1
  rsum=rsum+rdt(j)
enddo
rsumi=one/rsum
dsumi=rsumi/dble(nt)

 !azonal part:
do m=2,nt
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*d2lon(1,m)
  d0=-pm*tlat(1)

  deti=one/(a0*d0-b0*c0)
  ala0(1,m)=a0*deti
  bla0(1,m)=b0*deti
  cla0(1,m)=c0*deti
  dla0(1,m)=d0*deti

  cpla=pc12-0.1d0*d2lon(2,m)
  dpla=-0.1d0*tlat(2)
  ela1(1,m)= cpla*bla0(1,m)- apla*dla0(1,m)
  ela2(1,m)= dpla*bla0(1,m)-0.2d0*dla0(1,m)
  fla1(1,m)= apla*cla0(1,m)- cpla*ala0(1,m)
  fla2(1,m)=0.2d0*cla0(1,m)- dpla*ala0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-d2lon(j,m)
    d0=-tlat(j)
    cmla(j,m)=pc12-0.1d0*d2lon(j-1,m)
    dmla(j)=-0.1d0*tlat(j-1)

    ala0(j,m)=   ela1(j-1,m)*amla     +fla1(j-1,m)*0.2d0
    bla0(j,m)=b0+ela2(j-1,m)*amla     +fla2(j-1,m)*0.2d0
    cla0(j,m)=c0+ela1(j-1,m)*cmla(j,m)+fla1(j-1,m)*dmla(j)
    dla0(j,m)=d0+ela2(j-1,m)*cmla(j,m)+fla2(j-1,m)*dmla(j)
    deti=one/(ala0(j,m)*dla0(j,m)-bla0(j,m)*cla0(j,m))
    ala0(j,m)=ala0(j,m)*deti
    bla0(j,m)=bla0(j,m)*deti
    cla0(j,m)=cla0(j,m)*deti
    dla0(j,m)=dla0(j,m)*deti

    cpla=pc12-0.1d0*d2lon(j+1,m)
    dpla=-0.1d0*tlat(j+1)
    ela1(j,m)= cpla*bla0(j,m)- apla*dla0(j,m)
    ela2(j,m)= dpla*bla0(j,m)-0.2d0*dla0(j,m)
    fla1(j,m)= apla*cla0(j,m)- cpla*ala0(j,m)
    fla2(j,m)=0.2d0*cla0(j,m)- dpla*ala0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*d2lon(ng,m)
  d0=-pm*tlat(ng)

  cmla(ng,m)=pc12-0.1d0*d2lon(ngm1,m)
  dmla(ng)=-0.1d0*tlat(ngm1)

  ala0(ng,m)=a0+ela1(ngm1,m)*amla      +fla1(ngm1,m)*0.2d0
  bla0(ng,m)=b0+ela2(ngm1,m)*amla      +fla2(ngm1,m)*0.2d0
  cla0(ng,m)=c0+ela1(ngm1,m)*cmla(ng,m)+fla1(ngm1,m)*dmla(ng)
  dla0(ng,m)=d0+ela2(ngm1,m)*cmla(ng,m)+fla2(ngm1,m)*dmla(ng)

  deti=one/(ala0(ng,m)*dla0(ng,m)-bla0(ng,m)*cla0(ng,m))
  ala0(ng,m)=ala0(ng,m)*deti
  bla0(ng,m)=bla0(ng,m)*deti
  cla0(ng,m)=cla0(ng,m)*deti
  dla0(ng,m)=dla0(ng,m)*deti
enddo
!---------------------------------------------------------------------

return 
end subroutine

!===================================================================

subroutine laplinv(src,sol,der)
! Inverts Laplace's operator on src to give sol and tau times its
! latitudinal derivative der.  That is Lap(sol) = src
! and der = tau*d(sol)/dlat.
! *** All fields are in spectral space ***

implicit none

 !Passed variables:
double precision:: src(ng,nt),sol(ng,nt),der(ng,nt)
 !Local variables:
double precision:: rhs(ng),pfg(ng)
double precision:: gsum,const,pm,utdb,vtdb
integer:: j,m

!-------------------------------------------------------
 !Solve for the zonal part of der:
do j=1,ng
  rhs(j)=src(j,1)*rdt(j)
enddo
 !Above, rdt = rho/tau = dmu/dlat

gsum=f1112*(rhs(1)+rhs(ng))
do j=2,ngm1
  gsum=gsum+rhs(j)
enddo
const=gsum*rsumi

do j=1,ng
  rhs(j)=rhs(j)-const*rdt(j)
enddo

pfg(1)=dl24*(21.d0*rhs(1)+rhs(2))
do j=2,ngm1
  pfg(j)=pfg(j-1)+dl24*(rhs(j-1)+22.d0*rhs(j)+rhs(j+1))
enddo

rhs(1)=pfg(1)
do j=2,ngm1
  rhs(j)=pfg(j)+pfg(j-1)
enddo
rhs(ng)=pfg(ngm1)

der(1,1)=rhs(1)*htd0(1)

do j=2,ng
  der(j,1)=(rhs(j)-f14*der(j-1,1))*htd0(j)
enddo

do j=ngm1,1,-1
  der(j,1)=etd0(j)*der(j+1,1)+der(j,1)
enddo

do j=1,ng
  der(j,1)=der(j,1)*rhoi(j)
enddo

!-------------------------------------------------------
 !Solve for the zonal part of sol:
pfg(1)=dl24*(21.d0*der(1,1)+der(2,1))
do j=2,ngm1
  pfg(j)=pfg(j-1)+dl24*(der(j-1,1)+22.d0*der(j,1)+der(j+1,1))
enddo
pfg(ng)=pfg(ngm1)+dl24*(der(ngm1,1)+21.d0*der(ng,1))

rhs(1)=pfg(1)
do j=2,ng
  rhs(j)=pfg(j)+pfg(j-1)
enddo

sol(1,1)=rhs(1)*htd0(1)

do j=2,ng
  sol(j,1)=(rhs(j)-f14*sol(j-1,1))*htd0(j)
enddo

do j=ngm1,1,-1
  sol(j,1)=etd0(j)*sol(j+1,1)+sol(j,1)
enddo

 !Remove global mean value of sol:
gsum=f1112*(sol(1,1)*rdt(1)+sol(ng,1)*rdt(ng))
do j=2,ngm1
  gsum=gsum+sol(j,1)*rdt(j)
enddo
const=gsum*rsumi

do j=1,ng
  sol(j,1)=sol(j,1)-const
enddo

!-----------------------------------------------------------------
 !Loop over azonal longitudinal wavenumbers and solve the 2x2 
 !block tridiagonal problem:
do m=2,nt
   !Multiply source by 1/tau^2:
  do j=1,ng
    pfg(j)=src(j,m)*tsqi(j)
  enddo

  pm=psig(ksig(m))

  rhs(1)=pm*pfg(1)+0.1d0*pfg(2)
  do j=2,ngm1
    rhs(j)=pfg(j)+0.1d0*(pfg(j-1)+pfg(j+1))
  enddo
  rhs(ng)=0.1d0*pfg(ngm1)+pm*pfg(ng)

  sol(1,m)=-rhs(1)*bla0(1,m)
  der(1,m)= rhs(1)*ala0(1,m)

  do j=2,ng
    utdb=      -sol(j-1,m)*amla     -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cmla(j,m)-der(j-1,m)*dmla(j)
    sol(j,m)=utdb*dla0(j,m)-vtdb*bla0(j,m)
    der(j,m)=vtdb*ala0(j,m)-utdb*cla0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=ela1(j,m)*sol(j+1,m)+ela2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=fla1(j,m)*sol(j+1,m)+fla2(j,m)*der(j+1,m)+der(j,m)
  enddo

   !Multiply derivative by tau so that tau*d(sol)/dlat is returned:
  do j=1,ng
    der(j,m)=der(j,m)*tau(j)
  enddo

enddo

return
end subroutine

!==========================================================================

 !Main end module
end module     
