program zonal
!  -------------------------------------------------------------------------
!  |   Computes the zonal mean absolute vorticity, <zeta+f>, and the zonal |
!  |   mean zonal velocity <u> from data in zz.r4 (generated by caps/casl).|
!  |                                                                       |
!  |   Output is to the unformatted, direct-access file "avg.r4" which     |
!  |   contains single-precision numerical data of the form                |
!  |                        t, <zeta+f>, <u>                               |
!  |   with each record corresponding to a given time read from zz.r4.     |
!  |                                                                       |
!  |   The file ext-zonal.asc and max-zonal.asc contains the abs max       |
!  |   values of <zeta+f> and <u>.                                         |
!  -------------------------------------------------------------------------

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use spectral

real:: zz(ng,nt),zzz(ng)
real:: uu(ng,nt),zuu(ng)
real:: tr4,aspsqm1,rlat,zzzmax,zuumax
double precision:: wka(ng,nt),wkb(ng,nt)
double precision:: qq(ng,nt)

!---------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open input data file:
open(43,file='zz.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

 !Open output files:
open(22,file='avg.r4',form='unformatted', access='direct', &
                 &  status='replace',recl=4+8*ng)
open(51,file='ext-zonal.asc',status='replace')

!---------------------------------------------------------------
 !Define Coriolis frequency, f = -2*Omega*rho'*tau:
aspsqm1=asp**2-one
do j=1,ng
  rlat=(dble(j)-f12)*dl-hpi
  cof(j)=fpole*sin(rlat)/sqrt(one+aspsqm1*cos(rlat)**2)
enddo

zfac=one/dble(nt)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(43,rec=loop,iostat=iread) tr4,zz
  if (iread .ne. 0) exit 
  write(*,'(a,f9.2)') ' Processing t = ',tr4

  do i=1,nt
    do j=1,ng
      qq(j,i)=dble(zz(j,i))
    enddo
  enddo

   !Invert Laplace's operator on the PV anomaly:
  call laplinv(qq,wka,wkb)
   !Here the streamfunction psi is wka while wkb = tau*d(psi)/dlat.

   !Get physical space velocity:
  call revfft(ng,nt,wkb,trig,factors)

   !Copy into uu:
  do i=1,nt
    do j=1,ng
      uu(j,i)=-real(wkb(j,i))
    enddo
  enddo

   !Compute zonal averages:
  do j=1,ng
    zzz(j)=zero
    zuu(j)=zero
  enddo

  do i=1,nt
    do j=1,ng
      zzz(j)=zzz(j)+zz(j,i)
      zuu(j)=zuu(j)+uu(j,i)
    enddo
  enddo

  do j=1,ng
    zzz(j)=zzz(j)*zfac+cof(j)
    zuu(j)=zuu(j)*zfac
  enddo

   !Write diagnostic data for this time:
  write(22,rec=loop) tr4,zzz,zuu

   !Find max abs values and write to a file:
  zzzmax=abs(zzz(1))
  zuumax=abs(zuu(1))
  do j=2,ng
    zzzmax=max(zzzmax,abs(zzz(j)))
    zuumax=max(zuumax,abs(zuu(j)))
  enddo
  write(51,'(f9.2,1x,f15.8)') tr4,zzzmax,zuumax

enddo

 !Close all files:
close(43)
close(22)
close(51)

write(*,*)
write(*,*) ' The results are ready in avg.r4; to view the results type'
write(*,*) ' zonalview'
write(*,*)
write(*,*) ' The extremal values vs time are listed in ext-zonal.asc'

end program