!#########################################################################
!  Computes the imbalanced spectra of height, vorticity, divergence and 
!  acceleration divergence and writes the data to ispectra.asc.
!  The output can be viewed using ispec_view.

!  Also saves the balanced spectra, which can be viewed with bspec_view.

!  *** Must run dgbal first ***

!           Written 6/4/2018 by D G Dritschel @ St Andrews
!#########################################################################

program ispectra

 !Import contants, parameters and common arrays:
use constants
use spectral

implicit none

 !Physical arrays:
double precision::  hh(ng,ng), zz(ng,ng), dd(ng,ng), gg(ng,ng)
double precision:: bhh(ng,ng),bzz(ng,ng),bdd(ng,ng),bgg(ng,ng)
double precision::  aa(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng)

 !Other local variables:
double precision:: spec(0:ng),rms,rmsb,rmsi,t
integer:: loop,iread,k

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Open input data files:
open(32,file='evolution/dd.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(33,file='evolution/gg.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(34,file='evolution/hh.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(35,file='evolution/zz.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(42,file='evolution/bdd.r8',form='unformatted',access='direct', &
                              status='old',recl=nbytes)
open(43,file='evolution/bgg.r8',form='unformatted',access='direct', &
                              status='old',recl=nbytes)
open(44,file='evolution/bhh.r8',form='unformatted',access='direct', &
                              status='old',recl=nbytes)
open(45,file='evolution/bzz.r8',form='unformatted',access='direct', &
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
  read(32,rec=loop,iostat=iread) t,dd
  if (iread .ne. 0) exit 
  read(33,rec=loop) t,gg
  read(34,rec=loop) t,hh
  read(35,rec=loop) t,zz

   !Read balanced fields:
  read(42,rec=loop) t,bdd
  read(43,rec=loop) t,bgg
  read(44,rec=loop) t,bhh
  read(45,rec=loop) t,bzz

  write(*,'(a,f9.2)') ' *** Processing t = ',t

   !Compute 1d spectra of balanced fields:
  aa=bzz
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(51,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(51,'(2(1x,f12.8))') alk(k),spec(k)
  enddo

  aa=bdd
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(52,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(52,'(2(1x,f12.8))') alk(k),spec(k)
  enddo

  aa=bgg
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(53,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(53,'(2(1x,f12.8))') alk(k),spec(k)
  enddo

  aa=bhh
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(54,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(54,'(2(1x,f12.8))') alk(k),spec(k)
  enddo

   !Compute 1d spectra of imbalanced fields as well as various norms:
  aa=zz-bzz
  rms=sqrt(dsumi*sum(zz**2))
  rmsb=sqrt(dsumi*sum(bzz**2))
  rmsi=sqrt(dsumi*sum(aa**2))
  write(71,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(61,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(61,'(2(1x,f12.8))') alk(k),spec(k)
  enddo

  aa=dd-bdd
  rms=sqrt(dsumi*sum(dd**2))
  rmsb=sqrt(dsumi*sum(bdd**2))
  rmsi=sqrt(dsumi*sum(aa**2))
  write(72,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(62,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(62,'(2(1x,f12.8))') alk(k),spec(k)
  enddo

  aa=gg-bgg
  rms=sqrt(dsumi*sum(gg**2))
  rmsb=sqrt(dsumi*sum(bgg**2))
  rmsi=sqrt(dsumi*sum(aa**2))
  write(73,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(63,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(63,'(2(1x,f12.8))') alk(k),spec(k)
  enddo

  aa=hh-bhh
  rms=sqrt(dsumi*sum(hh**2))
  rmsb=sqrt(dsumi*sum(bhh**2))
  rmsi=sqrt(dsumi*sum(aa**2))
  write(74,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(64,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(64,'(2(1x,f12.8))') alk(k),spec(k)
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
write(*,*) ' Use spec_view to image the various spectra.'
write(*,*)
write(*,*) ' Norms of various fields (full, balanced and imbalanced)'
write(*,*) ' are in znorms.asc, dnorms.asc, gnorms.asc & hnorms.asc in the'
write(*,*) ' evolution subdirectory.'
write(*,*)

 !End main program
end program ispectra
!=======================================================================
