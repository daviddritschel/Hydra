program ranpv

!=====================================================
!         Sets up a random PV anomaly field
!=====================================================

 !Import necessary modules:
use constants
use spectral
use force

 !Declarations:
implicit none

double precision:: pslon(ng,nt),gslon(ng,nt)
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt),wkd(ng,nt)
double precision:: var(ng,nt)
double precision:: pslat(nt),gslat(nt)
double precision:: qbar(ng)
double precision:: eps,rk0i,fac,vrms,avrms,r

integer, dimension(:), allocatable :: seed
integer:: i,j,k,k0,jseed,ic,m

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

write(*,*) ' The PV is set equal to f + a random field with variance'
write(*,*) ' spectrum proportional to k^5*exp(-2k^2/k_0^2), with zero'
write(*,*) ' global average and an rms value equal to vrms.'
write(*,*)

write(*,*) ' Enter the rms PV anomaly relative to f_pole:'
read(*,*) eps
vrms=eps*fpole

write(*,*) ' Enter k_0:'
read(*,*) k0

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
rk0i=one/dble(k0)
do m=1,nt
  do j=1,ng
    pslon(j,m)=(rk0i*wave(m)*clati(j))**2
    gslon(j,m)=exp(-pslon(j,m))
  enddo
enddo
do k=1,nt
  pslat(k)=(rk0i*wave(k))**2
  gslat(k)=exp(-pslat(k))
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
    wka(i,k)=wka(i,k)*gslat(k)
    wkb(i,k)=wka(i,k)*pslat(k)
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
var=gslon*(pslon*wkc+wkd)

 !Return var to physical space:
call revfft(ng,nt,var,trig,factors)

 !Remove global mean:
call zeroavg(var)

 !Find rms value:
call getrms(var,avrms)

 !Normalise:
fac=vrms/avrms
var=fac*var

qbar=zero
do i=1,nt
  do j=1,ng
    qbar(j)=qbar(j)+var(j,i)**2
  enddo
enddo
qbar=sqrt(clati*qbar/dble(nt))
open(77,file='qbar_init.asc',status='replace')
do j=1,ng
  write(77,*) qbar(j),(dble(j)-f12)*dl-hpi
enddo
close(77)

 !Add f:
do i=1,nt
  do j=1,ng
    var(j,i)=var(j,i)+cof(j)
  enddo
enddo

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,var
close(20)

end program ranpv
