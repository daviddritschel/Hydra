program area
! ------------------------------------------------------------------------- !
! Computes the area in each contour level l, determined from ind(j).        !
!                                                                           !
! Writes area.asc listing t vs area(l), l = 1, ..., l_max.                  !
!                                                                           !
! Uses analytical form of the surface function mu(z), given in the paper    !
! Dritschel & Boatto, Proc. R. Soc. A 20140890                              !
!                                                                           !
! Written 15 February 2015 by D G Dritschel @ St Andrews                    !
! ------------------------------------------------------------------------- !

use parameters

implicit none

integer,parameter:: nlevm=2000
!nlevm: Max number of distinct contour levels

integer:: loc(-nlevm:nlevm),nlev

double precision,parameter:: zero=0.d0,  one=1.d0, two=2.d0, three=3.d0
double precision,parameter:: f12=one/two, f13=one/three, f16=one/6.d0
double precision,parameter:: pi=3.14159265358979323846264338327950d0
double precision,parameter:: twopi=two*pi
 ! p1,p2: points for 2-point Gaussian Quadrature:
double precision,parameter:: p1=0.21132486541d0,p2=0.78867513459d0
double precision,parameter:: tp1=two*p1,ttp1=three*p1/two
double precision,parameter:: tp2=two*p2,ttp2=three*p2/two

double precision,parameter:: bsq=asp**2

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
! Calculates the diagnostic, here area, for all times in the data

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
   !Compute area:
  call comparea(n,npt,loop,t)
  loop=loop+1
enddo

return
end subroutine

!==========================================================================

subroutine initialise
! Initialises constants and arrays & opens files for reading/writing

implicit none

double precision:: t
integer:: n,npt

!-----------------------------------------------------
 !Open input file:
open(80,file='cont/synopsis.asc',status='old')

 !Read first record to determine contour levels below:
read(80,*) t,n,npt
rewind 80

 !Get number of distinct contour levels (nlev):
call getlev(n,npt)

 !Open output file:
open(20,file='area.asc',status='replace')

return
end subroutine

!==========================================================================

subroutine getlev(n,npt)
! Initialises distinct contour levels

implicit none

integer:: np(n),i1(n),ind(n),lev(n)
integer:: n,npt,j,l,m
logical:: newlev

!-----------------------------------------------------
open(90,file='cont/index000',form='unformatted',status='old')
read(90) np,i1,ind
close(90)

nlev=1
loc(ind(1))=1

if (n .eq. 1) return

lev(1)=ind(1)
do j=2,n
  newlev=.true.
  m=ind(j)
  do l=1,nlev
    if (m .eq. lev(l)) newlev=.false.
  enddo
  if (newlev) then
    nlev=nlev+1
    lev(j)=nlev
    loc(m)=nlev
  endif      
enddo

return
end subroutine

!==========================================================================

subroutine finalise
! Closes all files and finishes execution

close(80)
close(20)

write(*,*)
write(*,*) ' t vs area is listed in area.asc'
write(*,*)

return
end subroutine

!==========================================================================

subroutine comparea(n,npt,loop,t)
! Calculates the vorticity field at time frame "loop" from the contours

implicit none

 !Area in each level:
double precision:: area(2*nlevm+1)

 !Local variables:
double precision:: x(npt),y(npt),z(npt)
double precision:: dx(npt),dy(npt),dz(npt)
double precision:: dsa(npt),dsb(npt),d(npt),e(npt)
double precision:: a(npt),b(npt),c(npt)
double precision:: t,eta,del,deta,ddel,xtmp,ytmp,ztmp
double precision:: bdif,p,afac,ax,ay,az,sx,sy,sz
double precision:: x1,y1,z1(npt),dx1,dy1,dz1,dlon1(npt)
double precision:: x2,y2,z2(npt),dx2,dy2,dz2,dlon2(npt)
double precision:: bc,pref,q,cth,mu0,mu1,mu2,sarea
integer:: np(n),i1(n),i2(n),ind(n)
integer:: next(npt)
integer:: n,npt,loop
integer:: i,ia,ib,ibeg,iend,j,l
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
  i2(j)=iend
  do i=ibeg,iend-1
    next(i)=i+1
  enddo
  next(iend)=ibeg
enddo 

!----------------------------------------------------------------
 !Compute the cubic interpolation coefficients:
do i=1,npt
  ia=next(i)
  dx(i)=x(ia)-x(i)
  dy(i)=y(ia)-y(i)
  dz(i)=z(ia)-z(i)
enddo

do i=1,npt
  dsa(i)=dx(i)**2+dy(i)**2+dz(i)**2
  e(i)=sqrt(dsa(i))
enddo

do ib=1,npt
  i=next(ib)
  dsb(i)=dsa(ib)
  a(i)=-dx(ib)
  b(i)=-dy(ib)
  c(i)=-dz(ib)
enddo

do i=1,npt
  if (dx(i)*a(i)+dy(i)*b(i)+dz(i)*c(i) .gt. zero) then 
     !Set curvature to zero at corners:
    d(i)=zero
  else
    d(i)=(x(i)*(dy(i)*c(i)-b(i)*dz(i))+ &
          y(i)*(dz(i)*a(i)-c(i)*dx(i))+ &
          z(i)*(dx(i)*b(i)-a(i)*dy(i)))/ &
          sqrt((a(i)*dsa(i)-dx(i)*dsb(i))**2+ &
               (b(i)*dsa(i)-dy(i)*dsb(i))**2+ &
               (c(i)*dsa(i)-dz(i)*dsb(i))**2+1.d-36)
  endif
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npt
  ia=next(i)
  dsb(i)=e(i)*(d(ia)+d(i))
  bdif  =e(i)*(d(ia)-d(i))
  a(i)=f16*bdif-f12*dsb(i)
  b(i)=f12*(dsb(i)-bdif)
  c(i)=f13*bdif
enddo

!----------------------------------------------------------------
 !Initialise areas:
do l=1,nlev
  area(l)=zero
enddo

 !Loop over contour segments and acculate areas:
do i=1,npt
  xtmp=y(i)*dz(i)-z(i)*dy(i)
  ytmp=z(i)*dx(i)-x(i)*dz(i)
  ztmp=x(i)*dy(i)-y(i)*dx(i)
  afac=sqrt((dx(i)**2+dy(i)**2+dz(i)**2)/(xtmp**2+ytmp**2+ztmp**2))
  ax=afac*xtmp
  ay=afac*ytmp
  az=afac*ztmp
  sx=dy(i)*az-dz(i)*ay
  sy=dz(i)*ax-dx(i)*az
  sz=dx(i)*ay-dy(i)*ax
   !1st Gaussian point:
  eta=p1*(a(i)+p1*(b(i)+p1*c(i)))
  del=f12*p1*(one-p1)
  deta=a(i)+tp1*(b(i)+ttp1*c(i))
  ddel=f12-p1
  x1=x(i)+p*dx(i)+eta*ax+del*sx
  y1=y(i)+p*dy(i)+eta*ay+del*sy
  z1(i)=z(i)+p*dz(i)+eta*az+del*sz
  dx1=dx(i)+deta*ax+ddel*sx
  dy1=dy(i)+deta*ay+ddel*sy
!  dz1=dz(i)+deta*az+ddel*sz
  dlon1(i)=(x1*dy1-y1*dx1)/(x1**2+y1**2)
   !2nd Gaussian point:
  eta=p2*(a(i)+p2*(b(i)+p2*c(i)))
  del=f12*p2*(one-p2)
  deta=a(i)+tp2*(b(i)+ttp2*c(i))
  ddel=f12-p2
  x2=x(i)+p*dx(i)+eta*ax+del*sx
  y2=y(i)+p*dy(i)+eta*ay+del*sy
  z2(i)=z(i)+p*dz(i)+eta*az+del*sz
  dx2=dx(i)+deta*ax+ddel*sx
  dy2=dy(i)+deta*ay+ddel*sy
!  dz2=dz(i)+deta*az+ddel*sz
  dlon2(i)=(x2*dy2-y2*dx2)/(x2**2+y2**2)
enddo

 !Calculate area depending on aspect ratio of ellipsoid:
if (asp .lt. one) then
  bc=sqrt(abs(one-bsq))
  pref=bsq/bc
  mu0=f12*(one + pref*log((one+bc)/asp))
  sarea=4.d0*pi*mu0
  do j=1,n
    l=loc(ind(j))
    do i=i1(j),i2(j)
      cth=z1(i)
      q=sqrt(bsq+(bc*cth)**2)
      mu1=f12*(q*cth + pref*log((q+bc*cth)/asp))
      cth=z2(i)
      q=sqrt(bsq+(bc*cth)**2)
      mu2=f12*(q*cth + pref*log((q+bc*cth)/asp))
      area(l)=area(l)+(mu0-mu1)*dlon1(i)+(mu0-mu2)*dlon2(i)
    enddo
  enddo
else if (asp .gt. one) then
  bc=sqrt(abs(bsq-one))
  pref=bsq/bc
  mu0=f12*(one + pref*asin(bc/asp))
  sarea=4.d0*pi*mu0
  do j=1,n
    l=loc(ind(j))
    do i=i1(j),i2(j)
      cth=z1(i)
      q=sqrt(bsq-(bc*cth)**2)
      mu1=f12*(q*cth + pref*asin(bc*cth/asp))
      cth=z2(i)
      q=sqrt(bsq-(bc*cth)**2)
      mu2=f12*(q*cth + pref*asin(bc*cth/asp))
      area(l)=area(l)+(mu0-mu1)*dlon1(i)+(mu0-mu2)*dlon2(i)
    enddo
  enddo
else
  sarea=4.d0*pi
  do j=1,n
    l=loc(ind(j))
    do i=i1(j),i2(j)
      area(l)=area(l)+(one-z1(i))*dlon1(i)+(one-z2(i))*dlon2(i)
    enddo
  enddo
endif

 !Normalise areas:
do l=1,nlev
  area(l)=f12*area(l)
enddo

 !Ensure all areas are < 1/2 of the total surface area in magnitude:
do l=1,nlev
  if (area(l) .gt. f12*sarea) then
    area(l)=sarea-area(l)
  else if (area(l) .lt. -f12*sarea) then
    area(l)=sarea+area(l)
  endif
enddo

write(20,'(1x,f12.5,10(1x,f14.9))') t,(area(l),l=1,nlev)
write(* ,'(1x,f12.5,10(1x,f14.9))') t,(area(l),l=1,nlev)

return
end subroutine

!==========================================================================

 !Main end program
end program
