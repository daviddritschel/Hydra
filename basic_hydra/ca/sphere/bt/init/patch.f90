program patch

! Sets up any number of approximately circular vortex patches (which 
! are perfectly circular only on a sphere) to compare with the motion 
! of point vortices.

! Written 17 Jan 2015 by dgd @ St Andrews
! Generalised to multiple patches on 20 March 2015 (96% solar eclipse!).

use contours

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: qc(ng,nt)

!------------------------------------------------------------
 !Initialisation:
call init_contours
aspsqm1=asp**2-one
pif=pi/180.d0

 !Set vorticity jump into each vortex to be 4*pi:
vor=four*pi

write(*,*) ' Number of vortex patches?'
read(*,*) n

write(*,*)
npt=0
do j=1,n
  write(*,'(a,i2)') ' ===> Specify the characteristics of vortex ',j
  write(*,*)

  write(*,*) '  Latitude of vortex centre (degrees)?'
  read(*,*) rlatc
  rlatc=rlatc*pif
  clatc=cos(rlatc)
  slatc=sin(rlatc)
  tauc=one/sqrt(one+aspsqm1*clatc**2)

  write(*,*) ' Longitude of vortex centre (degrees)?'
  read(*,*) rlonc
  rlonc=rlonc*pif
  clonc=cos(rlonc)
  slonc=sin(rlonc)

  write(*,*) ' Radius of vortex (in degrees latitude)?'
  read(*,*) raddeg
  radrad=raddeg*pif
  rad=sin(radrad)
  vlc=cos(radrad)
  xc=vlc*clatc
  zc=vlc*slatc

   !Set up contour:
  np(j)=nint(twopi*rad/(amu*ell))
  dalp=twopi/dble(np(j))
  do i=1,np(j)
    alp=dalp*(dble(i)-f12)
    calp=cos(alp)
    salp=sin(alp)
    xx=xc-rad*tauc*slatc*salp
    yy=rad*calp
    zz=zc+rad*tauc*clatc*salp
    fac=one/sqrt(xx**2+yy**2+zz**2)
    ii=npt+i
    x(ii)=fac*(xx*clonc-yy*slonc)
    y(ii)=fac*(yy*clonc+xx*slonc)
    z(ii)=fac*zz
    next(ii)=ii+1
  enddo
  ind(j)=1
  i1(j)=npt+1
  i2(j)=npt+np(j)
  next(i2(j))=i1(j)
  npt=i2(j)
enddo

 !Obtain gridded values (using default vorticity jump dq in parameters.f90):
call con2grid(qc)

 !Correct for prescribed vorticity jump (vor) into the vortex:
vrat=vor/dq
do i=1,nt
  do j=1,ng
    qc(j,i)=vrat*qc(j,i)
  enddo
enddo

 !Find maximum value of qc to work out mean value of original field:
qcmax=qc(1,1)
do i=1,nt
  do j=1,ng
    qcmax=max(qcmax,qc(j,i))
  enddo
enddo
qcbar=vor-qcmax

 !Surface area:
area=dl**2/dsumci

 !Summed circulation of all vortices (before removing average):
circ=qcbar*area

write(*,'(a,f14.10)') '  Surface area/(4*pi) = ',area/(4.d0*pi)
write(*,'(a,f14.10)') '  Total Circulation = ',circ

 !Choose output (depending on algorithm choice):
write(*,*)
write(*,*) ' What would you like to do?'
write(*,*) ' (1) write initial gridded field for caps, or'
write(*,*) ' (2) write initial contour for casl?'
read(*,*) iopt

if (iopt .eq. 1) then
   !Write initial absolute vorticity field:
  open(20,file='qq_init.r8',form='unformatted', &
      & access='direct',status='new',recl=2*nbytes)
  write(20,rec=1) zero,qc
  close(20)

else
   !Write initial contour:
  open(80,file='cont/synopsis.asc',status='new')
  write(80,'(1x,f12.5,1x,i9,1x,i10)') zero,n,npt

  open(81,file='cont/index000',form='unformatted',status='new')
  write(81) np(1:n),i1(1:n),ind(1:n)
  close(81)

  open(82,file='cont/nodes000',form='unformatted',status='new')
  write(82) x(1:npt),y(1:npt),z(1:npt)
  close(82)
endif

end program
