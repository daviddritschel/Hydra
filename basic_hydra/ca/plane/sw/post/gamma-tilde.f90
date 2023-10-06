!#########################################################################
!  Computes the gamma-tilde = gamma + 2J(u,v) - 2*delta^2.
!  Writes gt.r8

!           Written 22/5/2018 by D G Dritschel @ St Andrews
!#########################################################################

program gtilde

use common
use evolution

implicit none

 !Physical arrays (dd is in evolution already):
double precision:: aa(ng,ng),gg(ng,ng)
double precision:: gt(ng,ng),bgt(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng)

 !Other local variables:
double precision:: spec(0:ng),rms,rmsb,rmsi
integer:: iopt,loop,iread,k

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

write(*,*) ' Choose one of the following two options: find gamma_tilde for'
write(*,*) '   (1) the full fields, or'
write(*,*) '   (2) all fields (full, balanced & imbalanced)'
write(*,*) ' Option?'
read(*,*) iopt

 !Open input & output data files:
open(32,file='evolution/dd.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(33,file='evolution/gg.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(34,file='evolution/zz.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(51,file='evolution/gt.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
 !For gamma_tilde spectra:
open(61,file='spectra/tspec.asc',status='replace')  

if (iopt .eq. 2) then
  open(42,file='evolution/bdd.r8',form='unformatted',access='direct', &
                                status='old',recl=nbytes)
  open(43,file='evolution/bgg.r8',form='unformatted',access='direct', &
                                status='old',recl=nbytes)
  open(44,file='evolution/bzz.r8',form='unformatted',access='direct', &
                                status='old',recl=nbytes)
  open(52,file='evolution/bgt.r8',form='unformatted',access='direct', &
                                status='replace',recl=nbytes)
  open(53,file='evolution/igt.r8',form='unformatted',access='direct', &
                                status='replace',recl=nbytes)
   !For balanced and imbalanced gamma_tilde spectra:
  open(62,file='spectra/tbspec.asc',status='replace')
  open(63,file='spectra/tispec.asc',status='replace')
   !For rms norms:
  open(71,file='evolution/tnorms.asc',status='replace')
endif
  
 !Open output data file:

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(32,rec=loop,iostat=iread) t,dd
  if (iread .ne. 0) exit 

  read(33,rec=loop) t,gg
  read(34,rec=loop) t,zz

  write(*,'(a,f9.2)') ' *** Processing t = ',t

  !Obtain velocity field from dd & zz:
  aa=dd
  call ptospc(ng,ng,aa,ds,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*ds
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)

  aa=zz
  call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*wke
  call xderiv(ng,ng,hrkx,wke,wkc)
  call yderiv(ng,ng,hrky,wke,wkd)

  wka=wka-wkd
  wkb=wkb+wkc
  call spctop(ng,ng,wka,uu,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkb,vv,xfactors,yfactors,xtrig,ytrig)

  !Form gamma-tilde:
  call jacob(uu,vv,aa)
  aa=two*(aa-dd**2)
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  wka=filt*wka
  call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)
  gt=gg+aa

  !Write data:
  write(51,rec=loop) t,gt

   !Compute 1d spectra of full fields:
  aa=gt
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(61,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(61,'(2(1x,f12.8))') alk(k),spec(k)
  enddo
  rms=sqrt(dsumi*sum(gt**2))

  if (iopt .eq. 2) then
    !Also find gamma_tilde for the balanced fields:
    read(42,rec=loop) t,dd
    read(43,rec=loop) t,gg
    read(44,rec=loop) t,zz

    !Obtain velocity field from dd & zz:
    aa=dd
    call ptospc(ng,ng,aa,ds,xfactors,yfactors,xtrig,ytrig)
    wke=rlap*ds
    call xderiv(ng,ng,hrkx,wke,wka)
    call yderiv(ng,ng,hrky,wke,wkb)

    aa=zz
    call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
    wke=rlap*wke
    call xderiv(ng,ng,hrkx,wke,wkc)
    call yderiv(ng,ng,hrky,wke,wkd)

    wka=wka-wkd
    wkb=wkb+wkc
    call spctop(ng,ng,wka,uu,xfactors,yfactors,xtrig,ytrig)
    call spctop(ng,ng,wkb,vv,xfactors,yfactors,xtrig,ytrig)

    !Form gamma-tilde (balanced):
    call jacob(uu,vv,aa)
    aa=two*(aa-dd**2)
    call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
    wka=filt*wka
    call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)
    bgt=gg+aa

    !Write data:
    write(52,rec=loop) t,bgt

    aa=bgt
    call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,spec)
    spec=log10(spmf*spec+1.d-32)
    write(62,'(f12.5,1x,i5)') t,kmaxred
    do k=1,kmaxred
      write(62,'(2(1x,f12.8))') alk(k),spec(k)
    enddo
    rmsb=sqrt(dsumi*sum(bgt**2))

    aa=gt-bgt
    write(53,rec=loop) t,aa
    rmsi=sqrt(dsumi*sum(aa**2))
    call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,spec)
    spec=log10(spmf*spec+1.d-32)
    write(63,'(f12.5,1x,i5)') t,kmaxred
    do k=1,kmaxred
      write(63,'(2(1x,f12.8))') alk(k),spec(k)
    enddo
    write(71,'(f9.2,3(1x,e14.7))') t,rms,rmsb,rmsi

  endif
enddo

 !Close files:
close(32)
close(33)
close(34)
if (iopt .eq. 2) then
  close(42)
  close(43)
  close(44)
  close(61)
  close(62)
  close(63)
  close(71)
endif
close(51)

write(*,*)
write(*,*) ' gamma-tilde written to gt.r8 in the evolution subdirectory'
if (iopt .eq. 2) then
  write(*,*) ' with the balanced & imbalanced parts in bgt.r8 & igt.r8'
endif
write(*,*)

 !End main program
end program
!=======================================================================
