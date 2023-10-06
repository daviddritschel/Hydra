!#########################################################################
!  Computes the balanced and imbalanced spectra of height, vorticity, 
!  divergence and acceleration divergence and writes the data to 
!  ispectra.asc.  The output can be viewed using the script sv.

!  *** Must run dgbal first ***

!           Written 11/2/2021 by D G Dritschel @ St Andrews
!#########################################################################

program ispectra

 !Import contants, parameters and common arrays:
use constants
use spectral

implicit none

 !Physical arrays:
double precision::  hh(ng,nt), zz(ng,nt), dd(ng,nt), gg(ng,nt)
double precision:: bhh(ng,nt),bzz(ng,nt),bdd(ng,nt),bgg(ng,nt)
double precision::  aa(ng,nt)

 !Spectral arrays:
double precision:: wka(ng,nt)

 !Other local variables:
double precision:: spec(ng),rms,rmsb,rmsi,t
real:: qqr4(ng,nt),tr4
integer:: loop,iread,m

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Open input data files:
open(32,file='evolution/dd.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(33,file='evolution/gg.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(34,file='evolution/hh.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(35,file='evolution/zz.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(42,file='evolution/bdd.r4',form='unformatted',access='direct', &
                              status='old',recl=nbytes)
open(43,file='evolution/bgg.r4',form='unformatted',access='direct', &
                              status='old',recl=nbytes)
open(44,file='evolution/bhh.r4',form='unformatted',access='direct', &
                              status='old',recl=nbytes)
open(45,file='evolution/bzz.r4',form='unformatted',access='direct', &
                              status='old',recl=nbytes)

 !Open output data files:
open(51,file='spectra/zbspec.asc',status='replace')
open(52,file='spectra/dbspec.asc',status='replace')
open(53,file='spectra/gbspec.asc',status='replace')
open(54,file='spectra/hbspec.asc',status='replace')

open(61,file='spectra/zispec.asc',status='replace')
open(62,file='spectra/dispec.asc',status='replace')
open(63,file='spectra/gispec.asc',status='replace')
open(64,file='spectra/hispec.asc',status='replace')

open(71,file='evolution/znorms.asc',status='replace')
open(72,file='evolution/dnorms.asc',status='replace')
open(73,file='evolution/gnorms.asc',status='replace')
open(74,file='evolution/hnorms.asc',status='replace')

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0

   !Read original fields:
  read(32,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit
  t=real(tr4)
  dd=real(qqr4)
  read(33,rec=loop) tr4,qqr4
  gg=real(qqr4)
  read(34,rec=loop) tr4,qqr4
  hh=real(qqr4)
  read(35,rec=loop) tr4,qqr4
  zz=real(qqr4)

   !Read balanced fields:
  read(42,rec=loop) tr4,qqr4
  bdd=real(qqr4)
  read(43,rec=loop) tr4,qqr4
  bgg=real(qqr4)
  read(44,rec=loop) tr4,qqr4
  bhh=real(qqr4)
  read(45,rec=loop) tr4,qqr4
  bzz=real(qqr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',t

   !Compute 1d longitudinal power spectra of the balanced fields:
  aa=bzz
  call forfft(ng,nt,aa,trig,factors)
  call longspec(aa,spec)
  spec=log10(spec+1.d-32)
  write(51,'(f13.6,1x,i5)') t,ng
  do m=1,ng
    write(51,'(2(1x,f12.8))') alm(m),spec(m)
  enddo

  aa=bdd
  call forfft(ng,nt,aa,trig,factors)
  call longspec(aa,spec)
  spec=log10(spec+1.d-32)
  write(52,'(f13.6,1x,i5)') t,ng
  do m=1,ng
    write(52,'(2(1x,f12.8))') alm(m),spec(m)
  enddo

  aa=bgg
  call forfft(ng,nt,aa,trig,factors)
  call longspec(aa,spec)
  spec=log10(spec+1.d-32)
  write(53,'(f13.6,1x,i5)') t,ng
  do m=1,ng
    write(53,'(2(1x,f12.8))') alm(m),spec(m)
  enddo

  aa=bhh
  call forfft(ng,nt,aa,trig,factors)
  call longspec(aa,spec)
  spec=log10(spec+1.d-32)
  write(54,'(f13.6,1x,i5)') t,ng
  do m=1,ng
    write(54,'(2(1x,f12.8))') alm(m),spec(m)
  enddo

   !Compute 1d longitudinal power spectra of the imbalanced fields:
  aa=zz-bzz
  call getrms(zz, rms)
  call getrms(bzz,rmsb)
  call getrms(aa, rmsi)
  write(71,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi
  call forfft(ng,nt,aa,trig,factors)
  call longspec(aa,spec)
  spec=log10(spec+1.d-32)
  write(61,'(f13.6,1x,i5)') t,ng
  do m=1,ng
    write(61,'(2(1x,f12.8))') alm(m),spec(m)
  enddo

  aa=dd-bdd
  call getrms(dd, rms)
  call getrms(bdd,rmsb)
  call getrms(aa, rmsi)
  write(72,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi
  call forfft(ng,nt,aa,trig,factors)
  call longspec(aa,spec)
  spec=log10(spec+1.d-32)
  write(62,'(f13.6,1x,i5)') t,ng
  do m=1,ng
    write(62,'(2(1x,f12.8))') alm(m),spec(m)
  enddo

  aa=gg-bgg
  call getrms(gg, rms)
  call getrms(bgg,rmsb)
  call getrms(aa, rmsi)
  write(73,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi
  call forfft(ng,nt,aa,trig,factors)
  call longspec(aa,spec)
  spec=log10(spec+1.d-32)
  write(63,'(f13.6,1x,i5)') t,ng
  do m=1,ng
    write(63,'(2(1x,f12.8))') alm(m),spec(m)
  enddo

  aa=hh-bhh
  call getrms(hh, rms)
  call getrms(bhh,rmsb)
  call getrms(aa, rmsi)
  write(74,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi
  call forfft(ng,nt,aa,trig,factors)
  call longspec(aa,spec)
  spec=log10(spec+1.d-32)
  write(64,'(f13.6,1x,i5)') t,ng
  do m=1,ng
    write(64,'(2(1x,f12.8))') alm(m),spec(m)
  enddo

enddo

 !Close files:
close(32)
close(33)
close(34)
close(35)

close(42)
close(43)
close(44)
close(45)

close(51)
close(52)
close(53)
close(54)

close(61)
close(62)
close(63)
close(64)

close(71)
close(72)
close(73)
close(74)

write(*,*)
write(*,*) ' Use sv to image the various spectra.'
write(*,*)
write(*,*) ' Norms of various fields (full, balanced and imbalanced)'
write(*,*) ' are in znorms.asc, dnorms.asc, gnorms.asc & hnorms.asc in the'
write(*,*) ' evolution subdirectory.'
write(*,*)

 !End main program
end program ispectra
!=======================================================================
