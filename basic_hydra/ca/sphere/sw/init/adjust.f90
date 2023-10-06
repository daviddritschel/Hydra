program adjust

!============================================================
!  Sets up a random height anomaly field in a state of rest
!============================================================

 !Import spectral module:
use spectral

 !Declarations:
implicit none

double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt),wkd(ng,nt)
double precision:: var(ng,nt)
double precision:: eps,rksri,fac,hmax,ahmax,r

integer, dimension(:), allocatable :: seed
integer:: i,j,k,ksr,jseed,ic,m

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

write(*,*) ' The height anomaly is set to a random field with a variance'
write(*,*) ' spectrum proportional to k^5*exp(-2k^2/ksr^2), with zero'
write(*,*) ' global average and a maximum absolute value equal to h_max.'
write(*,*)

write(*,*) ' Enter h_max:'
read(*,*) hmax

write(*,*) ' Enter ksr:'
read(*,*) ksr

write(*,*) ' Enter a random seed (integer):'
read(*,*) jseed

 !Initialize random number generator:
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=jseed
do i=1,jseed
  call random_seed(put=seed)
enddo

 !Initialize squared wavenumber arrays used below:
rksri=one/dble(ksr)
do m=1,nt
  do j=1,ng
    plon(j,m)=(rksri*wave(m)*clati(j))**2
    glon(j,m)=exp(-plon(j,m))
  enddo
enddo
do k=1,nt
  plat(k)=(rksri*wave(k))**2
  glat(k)=exp(-plat(k))
enddo

!------------------------------------------------------------------
 !Generate random values between -1 and +1 at all grid points:
call random_number(wkc)
wkc=two*wkc-one

 !Create great circles:
do i=1,ng
  ic=i+ng
  do j=1,ng
    wka(i,j)     =wkc(j,i)
    wka(i,ntp1-j)=wkc(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(ng,nt,wka,trig,factors)

 !Apply latitudinal spectrum:
do k=1,nt
  do i=1,ng
    wka(i,k)=wka(i,k)*glat(k)
    wkb(i,k)=wka(i,k)*plat(k)
  enddo
enddo

 !Inverse FFT in latitude:
call revfft(ng,nt,wka,trig,factors)
call revfft(ng,nt,wkb,trig,factors)

 !Unpack arrays:
do i=1,ng
  ic=i+ng
  do j=1,ng
    wkc(j,i) =wka(i,j)
    wkc(j,ic)=wka(i,ntp1-j)
    wkd(j,i) =wkb(i,j)
    wkd(j,ic)=wkb(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(ng,nt,wkc,trig,factors)
call forfft(ng,nt,wkd,trig,factors)

 !Apply longitudinal spectrum and define var:
var=glon*(plon*wkc+wkd)

 !Return var to physical space:
call revfft(ng,nt,var,trig,factors)

 !Remove global mean:
call zeroavg(var)

 !Normalise:
ahmax=maxval(abs(var))
fac=hmax/ahmax
var=fac*var

 !Write initial height anomaly field:
open(20,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,var
close(20)

 !Define PV for a rest state:
do i=1,nt
  do j=1,ng
    var(j,i)=fpole*slat(j)/(one+var(j,i))
  enddo
enddo

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,var
close(20)

 !Write initial divergence field:
var=zero
open(20,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,var
close(20)

end program adjust
