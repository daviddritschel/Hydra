program grid
!---------------------------------------------------------------------------
!       Creates a gridded vorticity field on the ultra-fine grid
!       Writes output to the fine subdirectory
!---------------------------------------------------------------------------

use contours

implicit none

integer,parameter:: mgu=16,ntu=mgu*nt, ngu=mgu*ng
 !mgu:       Ultra-fine conversion grid ratio
 !ntu,ngu:   number of grid boxes in the x & y directions

 !Contour -> Ultra-fine Grid arrays:
double precision:: qa(0:ngu+1,ntu+1),t
double precision:: clonu(ntu),slonu(ntu)
double precision:: dlu,dlui

 !For removing average of ultra-fine vorticity in congen.f90:
double precision:: rdtu(ngu-1),dsumui

real rqa(ngu,ntu)
character(len=3):: pind

!----------------------------------------------------------------
 !Initialise:
call init_grid

 !Convert contours to gridded values (qa)
call con2ugrid

 !Write data (single precision) to fine subdirectory:
open(44,file='fine/qq'//pind//'.r4',form='unformatted', &
      access='stream',status='replace')
do i=1,ntu
  do j=1,ngu
    rqa(j,i)=real(qa(j,i))
  enddo
enddo
write(44) real(t),rqa
close(44)

write(*,*)
write(*,*) ' Image by typing the command'
write(*,*)
write(*,'(a,i5,1x,i5)') ' dataview fine/qq'//pind//'.r4 -ndim ',nxu,nyu
write(*,*)


contains 

!=============================================================

subroutine init_grid
 !Initialises constants and arrays needed in this programme

implicit none
double precision:: rlonu,aspsqm1,rhou,rsumu
integer:: loop,i,ibeg,iend,j

!-------------------------------------------------------------------
 !Initialise contour module:
call init_contours

 !Constants used contour->ultra-fine grid conversion:
dlu =twopi/dble(ntu)
dlui=dble(ntu)/(twopi+small)

do i=1,ntu
  rlonu=dlu*dble(i-1)-pi
  clonu(i)=cos(rlonu)
  slonu(i)=sin(rlonu)
enddo

 !For removing average vorticity on the ultra-fine grid:
aspsqm1=asp**2-one
do j=1,ngu-1
  rhou=cos(dble(j)*dlu-hpi)
  rdtu(j)=rhou*sqrt(one+aspsqm1*rhou**2)
enddo
rsumu=zero
do j=1,ngu-1
  rsumu=rsumu+rdtu(j)
enddo
dsumui=one/(rsumu*dble(ntu))

!-------------------------------------------------------------------
 !Select period to process:
write(*,'(2(a,i3))') ' There are ',nperiod, &
                   & ' periods ranging from 0 to ',nperiod
write(*,*) ' Time period to obtain ultra-fine gridded vorticity?'
read(*,*) loop
write(pind(1:3),'(i3.3)') loop

write(*,'(a,f12.5)') ' *** Imaging t = ',dt*dble(nsteps*loop)

 !Read contour data:
open(40,file='cont/synopsis.asc',status='old')
do i=0,loop
  read(40,*) t,n,npt
enddo
close(40)

open(40,file='cont/index'//pind,form='unformatted',status='old')
read(40) np(1:n),i1(1:n),ind(1:n)
close(40)

open(40,file='cont/nodes'//pind,form='unformatted',status='old')
read(40) x(1:npt),y(1:npt),z(1:npt)
close(40)

 !Reconstruct next array:
do j=1,n
  ibeg=i1(j)
  iend=ibeg+np(j)-1
  do i=ibeg,iend-1
    next(i)=i+1
  enddo
  next(iend)=ibeg
enddo 

return
end subroutine

!=============================================================

subroutine con2ugrid
 !Converts contours to an ultra-find grid as qa:

implicit none
double precision:: cx(npt),cy(npt),cz(npt),sq(npt)
double precision:: sig,rlatc,p,qasp,qanp,aspsqm1,rlatu,avqa
integer:: ntc(npt),ilm1(npt)
integer:: i,j,k,ka,jump,ioff,ncr

 !Initialise crossing information:
do k=1,npt
  ilm1(k)=int(dlui*(pi+atan2(y(k),x(k))))
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
  ntc(k)=ntc(k)-ntu*((2*ntc(k))/ntu)
  if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
  if (abs(cz(k)) .gt. zero) then
    cx(k)=cx(k)/cz(k)
    cy(k)=cy(k)/cz(k)
  endif
enddo

!----------------------------------------------------------------------
 !Initialise PV jump array:
do i=1,ntu
  do j=0,ngu+1
    qa(j,i)=zero
  enddo
enddo

 !Determine crossing indices:
do k=1,npt
  if (ntc(k) .ne. 0) then
    jump=sign(1,ntc(k))
    ioff=ntu+ilm1(k)+(1+jump)/2
    ncr=0
    do while (ncr .ne. ntc(k))
      i=1+mod(ioff+ncr,ntu)
      rlatc=dlui*(hpi+atan(cx(k)*clonu(i)+cy(k)*slonu(i)))
      j=int(rlatc)+1
      p=rlatc-dble(j-1)
      qa(j,i)=  qa(j,i)+(one-p)*sq(k)
      qa(j+1,i)=qa(j+1,i)+    p*sq(k)
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get PV values, at half latitudes, by sweeping through latitudes:
do i=1,ntu
  do j=2,ngu
    qa(j,i)=qa(j,i)+qa(j-1,i)
  enddo
enddo
 !Here, qa(j,i) stands for the PV at latitude j-1/2,
 !from j = 1, ..., ngu.

 !Determine unique polar values:
qasp=zero
qanp=zero
do i=1,ntu
  qasp=qasp+qa(1  ,i)
  qanp=qanp+qa(ngu,i)
enddo
qasp=qasp/dble(ntu)
qanp=qanp/dble(ntu)

 !Average half-grid PV to full grid:
do i=1,ntu
  qa(0,i)=qasp
  do j=1,ngu-1
    qa(j,i)=f12*(qa(j,i)+qa(j+1,i))
  enddo
  qa(ngu,i)=qanp
enddo

 !Remove global average:
avqa=zero
do i=1,ntu
  do j=1,ngu-1
    avqa=avqa+qa(j,i)*rdtu(j)
  enddo
enddo
avqa=avqa*dsumui

 !Also remove f if omega .ne. 0:
if (rotate) then
  aspsqm1=asp**2-one
  do i=1,ntu
    do j=0,ngu
      rlatu=dlu*dble(j)-hpi
      qa(j,i)=qa(j,i)-avqa-fpole*sin(rlatu)/sqrt(one+aspsqm1*cos(rlatu)**2)
    enddo
  enddo
else
  do i=1,ntu
    do j=0,ngu
      qa(j,i)=qa(j,i)-avqa
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

 !End main program
end program
!=======================================================================
