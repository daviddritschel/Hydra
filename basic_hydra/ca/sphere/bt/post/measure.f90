program zonal
!  -------------------------------------------------------------------------
!  |   Computes the vorticity measure or area fraction of each vorticity   |
!  |   level at a selected time frame from data in grid/zz.r4              |
!  |                                                                       |
!  |   Output is to the formatted file "pdf.asc" which lists zeta vs       |
!  |   log10(area fraction).                                               |
!  -------------------------------------------------------------------------

use constants
implicit none

integer,parameter:: nbmax=1000
real:: zz(ng,nt),rdt(ng),pdf(nbmax)
real:: aspsqm1,rho,t,zzmin,zzmax,db,dbi,spdf
integer:: i,j,k,nb,period

!---------------------------------------------------------------
 !Define dmu/dtheta = rho*tau:
aspsqm1=asp**2-one
do j=1,ng
  rho=cos((dble(j)-f12)*dl-hpi)
  rdt(j)=rho*sqrt(one+aspsqm1*rho**2)
enddo

!---------------------------------------------------------------
 !Open input data file:
open(44,file='grid/zz.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

 !Select time frame:
write(*,*) ' Time period (choose 0 for the initial one)?'
read(*,*) period

 !Read data:
read(44,rec=period+1) t,zz
close(44)
write(*,'(a,f12.5)') ' t = ',t

 !Work out min/max values of field:
zzmax=zz(1,1)
zzmin=zz(1,1)
do i=1,nt
  do j=1,ng
    zzmax=max(zzmax,zz(j,i))
    zzmin=min(zzmin,zz(j,i))
  enddo
enddo
write(*,'(a,f12.6)') ' zeta_min = ',zzmin
write(*,'(a,f12.6)') ' zeta_max = ',zzmax

 !Initialise pdf:
write(*,*) ' How many bins do you wish to form the pdf?'
read(*,*) nb
db=(zzmax-zzmin)/float(nb-1)
dbi=1./db

do k=1,nb
  pdf(k)=0.
enddo

 !Form pdf:
do i=1,nt
  do j=1,ng
    k=nint((zz(j,i)-zzmin)*dbi)+1
    pdf(k)=pdf(k)+rdt(j)
  enddo
enddo

 !Normalise so the sum is 1:
spdf=0.
do k=1,nb
  spdf=spdf+pdf(k)
enddo
spdf=1./spdf

do k=1,nb
  pdf(k)=spdf*pdf(k)
enddo

 !Open output file and write data:
open(22,file='pdf.asc',status='replace')
do k=1,nb
  write(22,'(f12.6,1x,f10.6)') zzmin+db*(float(k)-0.5),log10(pdf(k))
enddo
close(22)

write(*,*)
write(*,*) ' zeta vs log10(area fraction) is ready in pdf.asc'

end program
