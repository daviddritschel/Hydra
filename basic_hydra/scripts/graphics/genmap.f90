program genmap
 !Generates the colourmap used by the hydra flow simulation package

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: nx=256*256
double precision,parameter:: dx=1.d0/dble(nx)
 !Scale factors controlling hyperbolas:
double precision,parameter:: sc1=1.d0/8.d0,sc2=sc1/2.d0
 !Maximum colour or saturation level (< 1):
double precision,parameter:: cm=sqrt(3.d0)/2.d0

double precision:: r(0:nx),g(0:nx),b(0:nx)

xi=0.25d0/sc1
a1=cm/sqrt(2.d0*xi+xi**2)

xi=0.125d0/sc2
a2=cm/sqrt(2.d0*xi+xi**2)

!-------------------------------
 !Red:
do i=0,nx
  x=dx*dble(i)
  if (x .lt. 0.25d0) then
    r(i)=cm
  else if (x .lt. 0.375d0) then
    xi=(x-0.375d0)/sc2
    r(i)=a2*sqrt(xi**2-2.d0*xi)
  else if (x .lt. 0.75d0) then
    r(i)=0.d0
  else
    xi=(x-0.75d0)/sc1
    r(i)=a1*sqrt(2.d0*xi+xi**2)
  endif
enddo

!-------------------------------
 !Green:
do i=0,nx
  x=dx*dble(i)
  if (x .lt. 0.25d0) then
    xi=x/sc1
    g(i)=a1*sqrt(2.d0*xi+xi**2)
  else if (x .lt. 0.5d0) then
    g(i)=cm
  else if (x .lt. 0.75d0) then
    xi=(x-0.75d0)/sc1
    g(i)=a1*sqrt(xi**2-2.d0*xi)
  else
    g(i)=0.d0
  endif
enddo

!-------------------------------
 !Blue:
do i=0,nx
  x=dx*dble(i)
  if (x .lt. 0.375d0) then
    b(i)=0.d0
  else if (x .lt. 0.5d0) then
    xi=(x-0.375d0)/sc2
    b(i)=a2*sqrt(2.d0*xi+xi**2)
  else
    b(i)=cm
  endif
enddo

!------------------------------------------------
 !Write colourmap:
open(11,file='colourmap',status='unknown')
 !Reserve i=0 for white background colour:
write(11,'(3(1x,f9.7))') 0.99999,0.99999,0.99999
do i=1,nx
  write(11,'(3(1x,f9.7))') r(nx-i),g(nx-i),b(nx-i)
enddo
close(11)

!------------------------------------------------------------------------
 !Write also a shortened colourmap of 256 entries for use by char2ps.f90:
open(11,file='short-colourmap',status='unknown')
 !Reserve ic=0 for white background colour:
write(11,'(3(1x,i3))') 255,255,255
do ic=1,255
  i=256*ic
  write(11,'(3(1x,i3))') nint(256*r(nx-i)),nint(256*g(nx-i)),nint(256*b(nx-i))
enddo
close(11)

end program
