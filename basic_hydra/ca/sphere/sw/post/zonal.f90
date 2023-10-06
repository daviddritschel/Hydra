program zonal
!  -------------------------------------------------------------------------
!  |   Computes the zonal mean zonal velocity u, <u>, height anomaly, <h>, |
!  |   eddy kinetic energy <h[(u-<u>)^2+v^2]>/2, zonal mean PV, <q>, and   |
!  |   PV flux <v(q-<q>)> as a function of latitude from data in hh.r4,    |
!  |   dd.r4 and zz.r4 (generated previously by caps).                     |
!  |                                                                       |
!  |   View the results using the script zv                                |
!  |                                                                       |
!  |   Written => 12.02.20.21 <= by D G Dritschel @ St Andrews             |
!  -------------------------------------------------------------------------

 ! Import constants and parameters:
use constants
 ! Import spectral module:
use spectral

double precision:: hh(ng,nt),dd(ng,nt),zz(ng,nt)
double precision:: uu(ng,nt),vv(ng,nt),qq(ng,nt)
double precision:: uds(ng,nt),vds(ng,nt)
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt)
double precision:: zuu(ng),zhh(ng),zek(ng),zqq(ng),zvq(ng),phi(ng),zwk(ng)
double precision:: zuumax,zhhmax,zekmax,zqqmax,t
double precision:: ekin,epot,etot,angm
double precision:: zfac
real:: qqr4(ng,nt),tr4
integer:: i,j,loop,iread,m

!---------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

zfac=one/dble(nt)
do j=1,ng
  phi(j)=(dble(j)-f12)*dl-hpi
enddo


!---------------------------------------------------------------
 !Open all input data files:
open(41,file='evolution/hh.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(42,file='evolution/dd.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(43,file='evolution/zz.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open all output data files:
open(51,file='evolution/zu.asc',status='replace')
open(52,file='evolution/zh.asc',status='replace')
open(53,file='evolution/zk.asc',status='replace')
open(54,file='evolution/zq.asc',status='replace')
open(55,file='evolution/zf.asc',status='replace')
open(15,file='evolution/zecomp.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(41,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit
  t=dble(tr4)
  hh=dble(qqr4)
  read(42,rec=loop,iostat=iread) tr4,qqr4
  dd=dble(qqr4)
  read(43,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)
  write(*,'(a,f9.2)') ' Processing t = ',t

   !Compute velocity field:
  wkc=dd
  call forfft(ng,nt,wkc,trig,factors)
  call laplinv(wkc,wka,vds)
  call deriv(ng,nt,rk,wka,uds) 
  wkc=zz
  call forfft(ng,nt,wkc,trig,factors)
  call laplinv(wkc,wka,wkb)
  call deriv(ng,nt,rk,wka,wkc)
  do m=1,nt
    uu(:,m)=clati*uds(:,m)-wkb(:,m)
    vv(:,m)=vds(:,m)+clati*wkc(:,m)
  enddo
  call revfft(ng,nt,uu,trig,factors) 
  call revfft(ng,nt,vv,trig,factors) 

   !Define PV:
  do i=1,nt
    qq(:,i)=(cof+zz(:,i))/(one+hh(:,i))
  enddo
  
   !Compute zonal averages:
  zuu=zero
  zhh=zero
  zek=zero
  zqq=zero
  zvq=zero

  do j=1,ng
    zuu(j)=zuu(j)+zfac*sum(uu(j,:))
    zhh(j)=zhh(j)+zfac*sum(hh(j,:))
    zqq(j)=zqq(j)+zfac*sum(qq(j,:))
  enddo

   !Compute energy components and angular momentum:
  zwk=clat*(one+zhh)*zuu**2
  ekin=twopi*rsumi*(f1112*(zwk(1)+zwk(ng))+sum(zwk(2:ngm1)))
   !Assume here that zbb = 0 (bb not saved to compute its average):
  zwk=clat*zhh**2
  epot=twopi*csq*rsumi*(f1112*(zwk(1)+zwk(ng))+sum(zwk(2:ngm1)))
  etot=ekin+epot
  zwk=clat*(clat*((one+zhh)*zuu+omega*clat*zhh))
  angm=fourpi*rsumi*(f1112*(zwk(1)+zwk(ng))+sum(zwk(2:ngm1)))

   !Write energies & angular momentum to zecomp.asc:
  write(15,'(f13.6,4(1x,f16.9))') t,ekin,epot,etot,angm

  do i=1,nt
    wka(:,i)=f12*(one+hh(:,i))*((uu(:,i)-zuu)**2+vv(:,i)**2)
    wkb(:,i)=vv(:,i)*(qq(:,i)-zqq)
  enddo
  
  do j=1,ng
    zek(j)=zek(j)+zfac*sum(wka(j,:))
    zvq(j)=zvq(j)+zfac*sum(wkb(j,:))
  enddo

   !Write diagnostic data for this time:
  write(51,'(f13.6,1x,i5)') t,ng
  write(52,'(f13.6,1x,i5)') t,ng
  write(53,'(f13.6,1x,i5)') t,ng
  write(54,'(f13.6,1x,i5)') t,ng
  write(55,'(f13.6,1x,i5)') t,ng
  do j=1,ng
    write(51,'(2(1x,f12.8))') zuu(j),phi(j)
    write(52,'(2(1x,f12.8))') zhh(j),phi(j)
    write(53,'(2(1x,f12.8))') zek(j),phi(j)
    write(54,'(2(1x,f12.8))') zqq(j),phi(j)
    write(55,'(2(1x,f12.8))') zvq(j),phi(j)
  enddo

enddo

 !Close all files:
close(41)
close(42)
close(43)
close(51)
close(52)
close(53)
close(54)
close(55)
close(15)

write(*,*)
write(*,*) ' To view the results type'
write(*,*)
write(*,*) ' zv'
write(*,*)

end program
