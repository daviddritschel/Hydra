program decay
!  ---------------------------------------------------------------------------
!  |   From the PV enstrophy spectra in spectra.asc, this routine computes   |
!  |   the total potential enstrophy, kinetic energy and full energy, then   |
!  |   3 characteristic wavenumbers obtained from the first moments of the   |
!  |   potential enstrophy, kinetic energy and full energy spectra.          |
!  |                                                                         |
!  |   Output is to the ascii files "qke.asc" which contains t vs the        |
!  |   potential enstrophy, kinetic energy and full energy, and "kqke.asc"   |
!  |   which contains the respective characteristic wavenumbers.             |
!  |                                                                         |
!  |   Written 8 April 2016 by D G Dritschel @ St Andrews                    |
!  ---------------------------------------------------------------------------

 !Import constants and parameters:
use constants
use spectral

implicit none

 !Declarations:
integer,parameter:: nmax=max(nx,ny)
integer:: k,loop,iread,maxwave

double precision,parameter:: smf=two*pi**2
double precision:: wave(nmax),ksq(nmax),kksq(nmax)
double precision:: qspec(nmax),kspec(nmax),espec(nmax)
double precision:: t,dum,sumqspec,qql2
double precision:: qtot,ktot,etot
double precision:: qwave,kwave,ewave

 !Define wavenumbers for use in sums below:
do k=1,nmax
  wave(k)=dble(k)
  ksq(k)=wave(k)**2
  kksq(k)=ksq(k)+kdsq
enddo

!---------------------------------------------------------------
 !Open file containing potential enstrophy spectra at all times:
open(51,file='spectra.asc',status='old')

 !Open output files:
open(21,file= 'qke.asc',status='replace')
open(31,file='kqke.asc',status='replace')

 !Loop over times available and process:
loop=0
do  
  loop=loop+1

  iread=0
  read(51,*,iostat=iread) t,sumqspec,qql2,maxwave
  if (iread .ne. 0) exit 

  write(*,'(a,f12.5)') ' Processing t = ',t

  do k=1,maxwave
    read(51,*) dum,qspec(k)
  enddo

   !Form spectra (include 2*pi^2 factor to give correct L2 norm):
  do k=1,maxwave
    qspec(k)=smf*10.d0**qspec(k)
    kspec(k)=ksq(k)*qspec(k)/kksq(k)**2
    espec(k)=qspec(k)/kksq(k)
  enddo

   !Compute diagnostics:
  qtot=zero
  ktot=zero
  etot=zero

  qwave=zero
  kwave=zero
  ewave=zero

  do k=1,maxwave
    qtot=qtot+qspec(k)
    qwave=qwave+qspec(k)*wave(k)
    ktot=ktot+kspec(k)
    kwave=kwave+kspec(k)*wave(k)
    etot=etot+espec(k)
    ewave=ewave+espec(k)*wave(k)
  enddo

  qwave=qwave/qtot
  kwave=kwave/ktot
  ewave=ewave/etot

   !Write data:
  write(21,'(f9.2,3(1x,e14.7))') t, qtot, ktot, etot
  write(31,'(f9.2,3(1x,f14.9))') t,qwave,kwave,ewave

enddo

close(21)
close(31)
close(51)

 !End main program
end program
!=======================================================================
