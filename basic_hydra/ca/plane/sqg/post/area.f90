program conread

use constants

implicit none

double precision,parameter:: p1=0.21132486541d0,p2=one-p1
 !p1,p2: points for 2-point Gaussian Quadrature
double precision,parameter:: dp1=two*p1,tp1=1.5d0*p1
double precision,parameter:: dp2=two*p2,tp2=1.5d0*p2
integer,parameter:: nlev=1000

double precision:: x(npm),y(npm),dx(npm),dy(npm)
double precision:: a(npm),b(npm),c(npm)
double precision:: u(npm),v(npm),e(npm)
double precision:: area(-nlev:nlev)
double precision:: t,dq,x1,x2,dx1,dx2,y1,y2,dy1,dy2
double precision:: xx,yy,eta,deta

integer:: np(nm),i1(nm),i2(nm),ind(nm),next(npm)
integer:: loop,iread,n,npt,i,j,ibeg,iend,ia,ib
integer:: levmin,levmax,lev

character(len=3):: pind

!-----------------------------------------------------------------
open(40,file='cont/qqsynopsis.asc',status='old')
loop=0
write(*,*) 'Index   n       npt        t'
do
   read(40,*,iostat=iread) n,npt,t,dq
   if (iread .ne. 0) exit
   loop=loop+1
   write(*,'(i3,2x,i6,2x,i8,2x,f11.6)') loop,n,npt,t
enddo

write(*,*)
write(*,*) ' Which index do you want to process?'
read(*,*) loop
write(pind(1:3),'(i3.3)') loop

rewind(40)
do i=1,loop
   read(40,*) n,npt,t,dq
enddo
close(40)

open(40,file='cont/qqindex'//pind,form='unformatted',status='old')
read(40) np(1:n),i1(1:n),ind(1:n)
close(40)

open(40,file='cont/qqnodes'//pind,form='unformatted',status='old')
read(40) x(1:npt),y(1:npt)
close(40)

 !Reconstruct nextq array:
do j=1,n
   ibeg=i1(j)
   iend=ibeg+np(j)-1
   i2(j)=iend
   do i=ibeg,iend-1
      next(i)=i+1
   enddo
   next(iend)=ibeg
enddo

write(*,*)
write(*,*) ' Processing data...'

 !Compute the contour increments:
do i=1,npt
   ia=next(i)
   xx=x(ia)-x(i)
   dx(i)=xx-ellx*dble(int(xx*hlxi))
   yy=y(ia)-y(i)
   dy(i)=yy-elly*dble(int(yy*hlyi))
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npt
   v(i)=dx(i)*dx(i)+dy(i)*dy(i)
   e(i)=sqrt(v(i))
enddo

do ib=1,npt
   i=next(ib)
   u(i)=v(ib)
   a(i)=-dx(ib)
   c(i)=-dy(ib)
enddo

do i=1,npt
   if (dx(i)*a(i)+dy(i)*c(i) > zero) then 
      !Set curvature to zero at corners:
      b(i)=zero
   else
      b(i)=(dx(i)*c(i)-a(i)*dy(i))/ &
           sqrt((a(i)*v(i)-dx(i)*u(i))**2+(c(i)*v(i)-dy(i)*u(i))**2+small3)
   endif
enddo

do i=1,npt
   ia=next(i)
   u(i)=e(i)*(b(ia)+b(i))
   c(i)=e(i)*(b(ia)-b(i))
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npt
   a(i)=f16*c(i)-f12*u(i)
   b(i)=f12*(u(i)-c(i))
   c(i)=f13*c(i)
enddo

 !Levels to process:  
levmin=minval(ind)
levmax=maxval(ind)
area(levmin:levmax)=zero

 !Accumulate areas:
do j=1,n
   lev=ind(j)
   do i=i1(j),i2(j)
      eta=p1*(a(i)+p1*(b(i)+p1*c(i)))
      deta=a(i)+dp1*(b(i)+tp1*c(i))
      x1=x(i)+p1*dx(i)-eta*dy(i)
      y1=y(i)+p1*dy(i)+eta*dx(i)
      dx1=dx(i)-deta*dy(i)
      dy1=dy(i)+deta*dx(i)
      eta=p2*(a(i)+p2*(b(i)+p2*c(i)))
      deta=a(i)+dp2*(b(i)+tp2*c(i))
      x2=x(i)+p2*dx(i)-eta*dy(i)
      y2=y(i)+p2*dy(i)+eta*dx(i)
      dx2=dx(i)-deta*dy(i)
      dy2=dy(i)+deta*dx(i)
      area(lev)=area(lev)+x1*dy1-y1*dx1+x2*dy2-y2*dx2
   enddo
enddo

area(levmin:levmax)=f14*area(levmin:levmax)

 !Write data:
open(33,file='area'//pind//'.asc',status='replace')
do lev=levmin,-1
   write(33,'(2(1x,f15.11))') dq*(dble(lev)+f12),area(lev)
enddo

do lev=1,levmax
   write(33,'(2(1x,f15.11))') dq*(dble(lev)-f12),area(lev)
enddo
close(33)

write(*,*)
write(*,*) ' All done.  q vs area is in area'//pind//'.asc'

end program conread
