!#########################################################################
!  Computes the integrated pressure P = p_n + 0.5*c^2*(1+h)^2, and
!  optionally its balanced part P_b = 0.5*c^2*(1+h_b)^2.

!  Writes pp.r8 (and optionally bpp.r8) and files for norms and spectra.

!           Written 26/4/2021 by D G Dritschel @ St Andrews
!#########################################################################

program pressure

 !Import contants, parameters and common arrays:
use constants
use spectral

implicit none

 !Physical arrays:
double precision:: aa(ng,ng),hh(ng,ng),pp(ng,ng),bpp(ng,ng)

 !Spectral array:
double precision:: wka(ng,ng)

double precision:: spec(0:ng),rms,rmsb,rmsi,t
integer:: iopt,loop,iread,k

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

write(*,*) ' Choose one of the following two options: find P for'
write(*,*) '   (1) the full fields, or'
write(*,*) '   (2) all fields (full, balanced & imbalanced)'
write(*,*) ' Option?'
read(*,*) iopt

 !Open input & output data files:
open(31,file='evolution/hh.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(36,file='evolution/pn.r8',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(51,file='evolution/pp.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(60,file='spectra/pnspec.asc',status='replace') !replaces old pspec.asc
open(61,file='spectra/pspec.asc',status='replace')  

if (iopt .eq. 2) then
  open(41,file='evolution/bhh.r8',form='unformatted',access='direct', &
                                status='old',recl=nbytes)
  open(52,file='evolution/bpp.r8',form='unformatted',access='direct', &
                                status='replace',recl=nbytes)
  open(62,file='spectra/pbspec.asc',status='replace')
  open(63,file='spectra/pispec.asc',status='replace')
  open(71,file='evolution/pnorms.asc',status='replace')
endif

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) t,hh
  if (iread .ne. 0) exit 

  read(36,rec=loop) t,pp

  write(*,'(a,f9.2)') ' *** Processing t = ',t

   !First obtain p_n spectrum (to replace old pspec.asc):
  aa=pp
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(60,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(60,'(2(1x,f12.8))') alk(k),spec(k)
  enddo
  
   !Define full integrated pressure P:
  aa=hh*(one+f12*hh)
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  wka=filt*wka
  call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)
  pp=pp+csq*aa

  !Write data:
  write(51,rec=loop) t,pp

   !Compute 1d spectra of full field:
  aa=pp
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,spec)
  spec=log10(spmf*spec+1.d-32)
  write(61,'(f12.5,1x,i5)') t,kmaxred
  do k=1,kmaxred
    write(61,'(2(1x,f12.8))') alk(k),spec(k)
  enddo
  rms=sqrt(dsumi*sum(pp**2))

  if (iopt .eq. 2) then
     !Also obtain balanced P:
    read(41,rec=loop) t,hh

     !Define balanced integrated pressure P:
    aa=hh*(one+f12*hh)
    call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
    wka=filt*wka
    call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)
    bpp=csq*aa

    !Write data:
    write(52,rec=loop) t,bpp

    aa=bpp
    call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,spec)
    spec=log10(spmf*spec+1.d-32)
    write(62,'(f12.5,1x,i5)') t,kmaxred
    do k=1,kmaxred
      write(62,'(2(1x,f12.8))') alk(k),spec(k)
    enddo
    rmsb=sqrt(dsumi*sum(bpp**2))

    aa=pp-bpp
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
close(31)
close(36)
close(51)
close(60)
close(61)
if (iopt .eq. 2) then
  close(41)
  close(52)
  close(62)
  close(63)
  close(71)
endif

write(*,*)
write(*,*) ' P is written to pp.r8 in the evolution subdirectory'
if (iopt .eq. 2) then
  write(*,*) ' with the balanced part in bpp.r8'
endif
write(*,*)

 !End main program
end program pressure
!=======================================================================
