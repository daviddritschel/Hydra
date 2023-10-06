program fspec
!  -------------------------------------------------------------------------
!  |   Computes the potential enstrophy spectra for a selected range of    |
!  |   datasets in the fine subdirectory, containing the PV on the         |
!  |   ultra-fine grid.                                                    |
!  |                                                                       |
!  |   Output is to the formatted files "fine/specxxx" where xxx is the    |
!  |   three digit identifier for the data set (the "period").             |
!  -------------------------------------------------------------------------

 !Import constants and parameters:
use constants
use sta2dfft

implicit none

 !Maximum x & y wavenumbers:
integer,parameter:: nwx=nxu/2,nwy=nyu/2

double precision:: qq(nyu,nxu),ss(nxu,nyu)
double precision::  alk(0:max(nxu,nyu))
double precision:: spmf(0:max(nxu,nyu)),spec(0:max(nxu,nyu))
double precision:: rkx(nxu),hrkx(nxu)
double precision:: rky(nyu),hrky(nyu)
double precision:: xtrig(2*nxu),ytrig(2*nyu)
double precision:: scx,rkxmax
double precision:: scy,rkymax
double precision:: delk,delki,snorm

real:: rqq(nyu,nxu),tr4

integer:: kmag(nxu,nyu),xfactors(5),yfactors(5)
integer:: kx,ky,k,kmax,kmaxred
integer:: loop1,loop2,loop
character(len=3):: pind

 !-----------------------------------------------------------------------
 !Set up FFT trig tables:
write(*,*) ' Initialising trig tables for FFTs...'
call init2dfft(nxu,nyu,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define x wavenumbers:
rkx(1)=zero
do kx=1,nwx-1
  rkx(kx+1)    =hrkx(2*kx)
  rkx(nxu+1-kx)=hrkx(2*kx)
enddo
rkx(nwx+1)=hrkx(nxu)
scx=twopi/ellx
rkxmax=scx*dble(nwx)

 !Define y wavenumbers:
rky(1)=zero
do ky=1,nwy-1
  rky(ky+1)    =hrky(2*ky)
  rky(nyu+1-ky)=hrky(2*ky)
enddo
rky(nwy+1)=hrky(nyu)
scy=twopi/elly
rkymax=scy*dble(nwy)

 !Initialise fixed arrays for computing spectra:
delk=sqrt(f12*(scx**2+scy**2))
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do ky=1,nyu
  do kx=1,nxu
    k=nint(sqrt(rkx(kx)**2+rky(ky)**2)*delki)
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(nxu*nyu) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi/dble(nxu*nyu)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(delk*dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=nint(sqrt((rkxmax**2+rkymax**2)/two)*delki)

 !-----------------------------------------------------------------------
 !Select datasets and process data:
write(*,*) ' Range of periods to process?'
read(*,*) loop1,loop2

do loop=loop1,loop2

  write(pind(1:3),'(i3.3)') loop

   !Open input file:
  open(44,file='fine/qq'//pind//'.r4',form='unformatted', &
        access='stream',status='old')

   !Read data:
  read(44) tr4,rqq
  close(44)
  write(*,'(a,f12.5)') '  Processing t = ',tr4

   !Convert to double precision:
  qq=dble(rqq)

   !Take a Fourier transform:
  call ptospc(nxu,nyu,qq,ss,xfactors,yfactors,xtrig,ytrig)

   !Calculate spectrum:
  do k=0,kmax
    spec(k)=zero
  enddo

   !x and y-independent mode:
  k=kmag(1,1)
  spec(k)=spec(k)+f14*ss(1,1)**2

   !y-independent mode:
  do kx=2,nxu
    k=kmag(kx,1)
    spec(k)=spec(k)+f12*ss(kx,1)**2
  enddo

   !x-independent mode:
  do ky=2,nyu
    k=kmag(1,ky)
    spec(k)=spec(k)+f12*ss(1,ky)**2
  enddo

   !All other modes:
  do ky=2,nyu
    do kx=2,nxu
      k=kmag(kx,ky)
      spec(k)=spec(k)+ss(kx,ky)**2
    enddo
  enddo

   !Normalise to take into account uneven sampling of wavenumbers 
   !in each shell [k-1/2,k+1/2]:
  do k=1,kmax
    spec(k)=spmf(k)*spec(k)
  enddo

   !Write data:
  open(51,file='fine/spec'//pind,status='replace')
  do k=1,kmaxred
    write(51,'(2(1x,f12.8))') alk(k),log10(spec(k))
  enddo
  close(51)

enddo

write(*,*)
write(*,*) ' All done.  See fine/specxxx for the PV spectra.'

 !End main program
end program
!=======================================================================
