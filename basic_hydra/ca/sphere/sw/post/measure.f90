program measure
!  -------------------------------------------------------------------------
!  |   Computes the equivalent latitude as a function of PV.  Output is    |
!  |   to the file evolution/qe.asc, in the same format used in zonal.f90  |
!  -------------------------------------------------------------------------

use constants
implicit none

real:: qq(ng,nt),dar(ng),area(0:ng)
real:: t,qqmin,qqmax,db,dbi
integer:: i,j,loop,iread,k

!---------------------------------------------------------------
!Define cell-centred area / 2*pi of each grid cell:
db=dl/twopi
do j=1,ng
  dar(j)=db*(sin(float(j)*dl-hpi)-sin(float(j-1)*dl-hpi))
enddo

!---------------------------------------------------------------
 !Open input data file:
open(31,file='evolution/qq.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output data file:
open(41,file='evolution/qe.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) t,qq
  if (iread .ne. 0) exit

  write(*,'(a,f9.2)') ' Processing t = ',t

   !Work out min/max values of field and divide into bins:
  qqmax=maxval(qq)
  qqmin=minval(qq)
  db=(qqmax-qqmin)/float(ng)
  dbi=1./db

   !Initialise area of each bin:
  area=0.

   !Accumulate areas:
  do i=1,nt
    do j=1,ng
      k=min(ng,int((qq(j,i)-qqmin)*dbi)+1)
      area(k)=area(k)+dar(j)
    enddo
  enddo

   !Output as equivalent latitude:
  write(41,'(f13.6,1x,i5)') t,ng
  write(41,'(2(1x,f12.8))') qqmin,-hpi
  area(0)=-1.
  do k=1,ng
    area(k)=area(k)+area(k-1)
    write(41,'(2(1x,f12.8))') qqmin+db*float(k),asin(area(k))
  enddo
enddo
close(31)
close(41)

write(*,*)
write(*,*) ' PV vs equivalent latitude is in evolution/qe.asc'

end program measure
