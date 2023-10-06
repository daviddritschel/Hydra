program enstro
! ------------------------------------------------------------------------- !
! Computes the enstrophy <q^2>/2 from contour data in the cont subdirectory !
!                                                                           !
! Writes ens.asc listing t vs enstrophy.                                    !
! ------------------------------------------------------------------------- !

use parameters

implicit none

integer,parameter:: mgf=16, ngf=ng*mgf, ntf=nt*mgf
!mgf:      Fine grid/inversion grid ratio used for the calculation here
!ntf,ngf:  Nnumber of grid boxes in the x & y directions on the fine grid

double precision,parameter:: zero=0.d0,  one=1.d0, two=2.d0
double precision,parameter:: f12=one/two, f1112=11.d0/12.d0
double precision,parameter:: pi=3.14159265358979323846264338327950d0
double precision,parameter:: hpi=pi/two,twopi=two*pi
double precision,parameter:: small=1.d-12

!Common fixed arrays and constants:
double precision:: clonf(ntf),slonf(ntf),rdtc(ngf)
double precision:: dlf,dlfi,dsumci

!---------------------------------------------------------------
 !Initialise:
call initialise

 !Compute enstrophy:
call diagnose

 !Finalise:
call finalise


!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains

!==========================================================================

subroutine diagnose
! Calculates the diagnostic, here enstrophy, for all times in the data

implicit none

double precision:: t
integer:: n,npt
integer:: loop,iread

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  iread=0
  read(80,*,iostat=iread) t,n,npt
  if (iread .ne. 0) exit 
   !Convert contours to gridded values:
  call con2grid(n,npt,loop,t)
  loop=loop+1
enddo

return
end subroutine

!==========================================================================

subroutine initialise
! Initialises constants and arrays & opens files for reading/writing

implicit none

double precision:: aspsqm1,rlonf,rhoc,rsumc,rsumci
integer:: i,j

!----------------------------------------------------------
 !Contour->grid conversion constants and arrays:
dlf =twopi/dble(ntf)
dlfi=dble(ntf)/(twopi+small)

do i=1,ntf
  rlonf=dlf*dble(i-1)-pi
  clonf(i)=cos(rlonf)
  slonf(i)=sin(rlonf)
enddo

!-----------------------------------------------------
 !For removing average vorticity:
aspsqm1=asp**2-one
do j=1,ngf
  rhoc=cos((dble(j)-f12)*dlf-hpi)
  rdtc(j)=rhoc*sqrt(one+aspsqm1*rhoc**2)
enddo

rsumc=f1112*(rdtc(1)+rdtc(ngf))
do j=2,ngf-1
  rsumc=rsumc+rdtc(j)
enddo
rsumci=one/rsumc
dsumci=rsumci/dble(ntf)

!-----------------------------------------------------
 !Open input file:
open(80,file='cont/synopsis.asc',status='old')

 !Open output file:
open(20,file='ens.asc',status='replace')

return
end subroutine

!==========================================================================

subroutine finalise
! Closes all files and finishes execution

close(80)
close(20)

write(*,*)
write(*,*) ' t vs enstrophy is listed in ens.asc'
write(*,*)

return
end subroutine

!==========================================================================

subroutine con2grid(n,npt,loop,t)
! Calculates the vorticity field at time frame "loop" from the contours

implicit none

!Fine grid vorticity field:
double precision:: zz(0:ngf+1,ntf)
 !Local variables:
double precision:: x(npt),y(npt),z(npt)
double precision:: cx(npt),cy(npt),cz(npt),sq(npt)
double precision:: t,sig,rlatc,p,vsum
integer:: np(n),i1(n),ind(n)
integer:: next(npt),ilm1(npt),ntc(npt)
integer:: n,npt,loop
integer:: i,ibeg,iend,j,k,ka
integer:: jump,ioff,ncr
character(len=3):: pind

!----------------------------------------------------------------
 !Open and read data files containing contour indices and nodes:
write(pind(1:3),'(i3.3)') loop
open(90,file='cont/index'//pind,form='unformatted',status='old')
read(90) np,i1,ind
close(90)

open(90,file='cont/nodes'//pind,form='unformatted',status='old')
read(90) x,y,z
close(90)

!----------------------------------------------------------------
 !Define next array for use below:
do j=1,n
  ibeg=i1(j)
  iend=ibeg+np(j)-1
  do i=ibeg,iend-1
    next(i)=i+1
  enddo
  next(iend)=ibeg
enddo 

!----------------------------------------------------------------
 !Initialise grid-line crossing information:
do k=1,npt
  ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
enddo

do k=1,npt
  ka=next(k)
  cx(k)=z(k)*y(ka)-y(k)*z(ka)
  cy(k)=x(k)*z(ka)-z(k)*x(ka)
  cz(k)=x(k)*y(ka)-y(k)*x(ka)
  ntc(k)=ilm1(ka)-ilm1(k)
enddo

do k=1,npt
  sig=sign(one,cz(k))
  sq(k)=dq*sig
  ntc(k)=ntc(k)-ntf*((2*ntc(k))/ntf)
  if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
  if (abs(cz(k)) .gt. zero) then
    cx(k)=cx(k)/cz(k)
    cy(k)=cy(k)/cz(k)
  endif
enddo

!----------------------------------------------------------------------
 !Initialise vorticity jump array:
do i=1,ntf
  do j=0,ngf+1
    zz(j,i)=zero
  enddo
enddo

 !Determine crossing indices:
do k=1,npt
  if (ntc(k) .ne. 0) then
    jump=sign(1,ntc(k))
    ioff=ntf+ilm1(k)+(1+jump)/2
    ncr=0
    do while (ncr .ne. ntc(k))
      i=1+mod(ioff+ncr,ntf)
      rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
      j=int(rlatc)+1
      p=rlatc-dble(j-1)
      zz(j,i)=  zz(j,i)+(one-p)*sq(k)
      zz(j+1,i)=zz(j+1,i)+    p*sq(k)
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get vorticity values, at half latitudes, by sweeping through latitudes:
do i=1,ntf
  do j=2,ngf
    zz(j,i)=zz(j,i)+zz(j-1,i)
  enddo
enddo
 !Here, zz(j,i) stands for the vorticity at latitude j-1/2,
 !from j = 1, ..., ngf.

!----------------------------------------------------------------
 !Remove average vorticity (compute 4th-order average):
vsum=zero
do i=1,ntf
  vsum=vsum+f1112*(rdtc(1)*zz(1,i)+rdtc(ngf)*zz(ngf,i))
  do j=2,ngf-1
    vsum=vsum+rdtc(j)*zz(j,i)
  enddo
enddo
vsum=vsum*dsumci

do i=1,ntf
  do j=1,ngf
    zz(j,i)=zz(j,i)-vsum
  enddo
enddo

!----------------------------------------------------------------
 !Compute enstrophy and write to a file:
vsum=zero
do i=1,ntf
  vsum=vsum+f1112*(rdtc(1)*zz(1,i)**2+rdtc(ngf)*zz(ngf,i)**2)
  do j=2,ngf-1
    vsum=vsum+rdtc(j)*zz(j,i)**2
  enddo
enddo
vsum=f12*vsum*dsumci

write(20,'(1x,f12.5,1x,f14.9)') t,vsum
write(* ,'(1x,f12.5,1x,f14.9)') t,vsum

return
end subroutine

!==========================================================================

 !Main end program
end program
