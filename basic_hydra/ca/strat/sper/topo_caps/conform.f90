program conform

use sta2dfft
use deriv1d

! ===================================================================
! This routine sets up the conformal transformation of the domain, 
! given the bottom topography h(x) in the function hp below.
! Modify that function as necessary.

! Use precomp to pass the parameter values for nx, ny, ellxp, ymin & ymax.
! Then compile using
! gfortran -O3 -o conform conform.f90 ~/hydra/lib/stafft/stafft.f90 ~/hydra/lib/stafft/sta2dfft.f90 ~/hydra/lib/stafft/deriv1d.f90 

! This is usually done through the "flow_setup" script.

! Adapted on 23 Aug 2015 by dgd from the python code topo_final.py 
! originally written by sek.
! ===================================================================

implicit double precision (a-h,o-z)

integer,parameter:: nx=N_X,ny=N_Y
integer,parameter:: nxm1=nx-1
integer,parameter:: nwx=nx/2,nwxm1=nwx-1
integer,parameter:: ngridp=nx*(ny+1),nbytes=4*(ngridp+1)
double precision,parameter:: ellxp=L_X,ymin=Y_MIN,ymax=Y_MAX
double precision,parameter:: elly=ymax-ymin
double precision,parameter:: tol=1.d-14
! tol: error tolerance in the conformal domain length

 !Local work arrays:
double precision:: x(0:nxm1),xp(0:nxm1),xpo(0:nxm1)
double precision:: gh(0:nxm1),h(0:nxm1),d(0:nxm1)
double precision:: rkx(0:nxm1),hrkx(nx),rky(ny)
double precision:: yori(0:ny,0:nxm1),xori(0:ny,0:nxm1)
double precision:: rkx1d(0:nxm1),dhdx(0:nxm1)
double precision:: yh0(0:ny)

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

!-------------------------------------------------------------------
! Initialise:
write(*,'(2(a,f6.2))') ' We consider a domain of height ',elly, &
                                           ' and length ',ellxp
write(*,*)
write(*,*) ' We take the topography to be of the form'
write(*,*) '       h(x) = A*exp(-b*(x-c)^2)'
write(*,*) ' [note: the domain is centred at x = 0.]'
write(*,*)
write(*,*) ' Enter A, b & c:'
read(*,*) a,b,c
! These are used in function hp below.

! Note, quantities relevant to the conformally mapped domain will NOT
! have a "p" appended to the end of their name.

! Initialise transformed space length (this is just a starting guess):
ellx=ellxp
hlx=ellx/2.d0
glx=ellx/dble(nx)

! Evenly-spaced conformal x coordinates at bottom of conformal domain:
do ix=0,nxm1
  x(ix)=glx*dble(ix)-hlx
enddo

! Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,rky)

! Define x wavenumbers (these are redefined below as ellx changes):
rkx(0)=0.d0
do kx=1,nwxm1
  rkx(kx)=hrkx(2*kx)
  rkx(nx-kx)=rkx(kx)
enddo
rkx(nwx)=hrkx(nx)

! Define spectral operator needed in inner loop below:
gh(0)=0.d0
do kx=1,nxm1
  gh(kx)=1.d0/tanh(rkx(kx)*elly)
enddo

! In outer iteration, adjust length (ellx) of conformal domain to 
! be consistent with the physical domain dimensions (ellxp,elly):
xerr=0.5d0
do while (abs(xerr) .gt. tol .and. abs(xerr) .lt. 1.d0)

  ! Ensure correct value for topography h(x) in the conformal domain:
  xperr=0.5d0
  xp=x
  do while (xperr .gt. tol .and. xperr .lt. 1.d0)

    !Assign h to hp(xp) and compute average (have):
    have=0.d0
    do ix=0,nxm1
      xpo(ix)=xp(ix)
      h(ix)=hp(xp(ix))
      have=have+h(ix)
    enddo
    have=have/dble(nx)

    ! FFT h:
    call forfft(1,nx,h,xtrig,xfactors)

    ! Multiply by i:
    d(0)=0.d0
    d(nwx)=0.d0
    do kx=1,nwxm1
      kxc=nx-kx
      d(kx)=-h(kxc)
      d(kxc)=h(kx)
    enddo

    ! Multiply by spectral gp operator:
    do kx=1,nxm1
      d(kx)=gh(kx)*d(kx)
    enddo

    ! Inverse FFT d:
    call revfft(1,nx,d,xtrig,xfactors)

    ! Complete new estimate of xp:
    xfac=1.d0-have/elly
    xp=xfac*x+d

    ! Compute error:
    terr=0.d0
    do ix=0,nxm1
      terr=terr+abs(xp(ix)-xpo(ix))
    enddo
    xperr=terr/dble(nx)

  enddo

  ! Correct ellx:
  xerr=(ellxp+2.d0*xp(0))/ellxp
  ellx=ellx+ellxp*xerr
  hlx=ellx/2.d0
  glx=ellx/dble(nx)

  ! Rescale wavenumbers:
  sfac=ellxp/ellx
  do kx=1,nwxm1
    rkx(kx)=sfac*hrkx(2*kx)
    rkx(nx-kx)=rkx(kx)
  enddo
  rkx(nwx)=sfac*hrkx(nx)

  ! Evenly-spaced conformal x coordinates at bottom of conformal domain:
  do ix=0,nxm1
    x(ix)=glx*dble(ix)-hlx
  enddo

  ! Re-define spectral operator needed in inner loop above:
  do kx=1,nxm1
    gh(kx)=1.d0/tanh(rkx(kx)*elly)
  enddo

enddo

if (abs(xerr) .le. tol) then

  ! Finalise values of h(x) and its average:
  have=0.d0
  do ix=0,nxm1
    h(ix)=hp(xp(ix))
    have=have+h(ix)
  enddo
  have=have/dble(nx)

  ! Write the converged value for h(x) to file:
  open(12,file='hh.asc',status='replace')
  do ix=0,nxm1
    write(12,'(1x,f17.14)') h(ix)
  enddo
  close(12)

  !------------------------------------------------------------
  ! Compute X and Y throughout the domain and write to a file:
  ! First set up 1D derivatives:
  call init_deriv(nx,ellx,rkx1d)

   !Form arrays for conformal mapping:

   !Fractional y grid values: 
  fac=1.d0/dble(ny)
  do iy=0,ny
    yh0(iy)=1.d0-fac*dble(iy)
  enddo
  gly=elly/dble(ny)

   !1D FFTs of h and dh/dx:
  call forfft(1,nx,h,xtrig,xfactors)
  call deriv(1,nx,rkx1d,h,dhdx)

   !Define the interior semi-spectral fields of X and Y:
  do iy=0,ny
    xori(iy,0)=0.d0
    yori(iy,0)=h(0)*yh0(iy)
  enddo
  do kx=1,nxm1
    fac=rkx(kx)*elly
    divy=1.d0/sinh(fac)
    divx=divy/rkx(kx)
    do iy=0,ny
      xori(iy,kx)=dhdx(kx)*cosh(fac*yh0(iy))*divx
      yori(iy,kx)=   h(kx)*sinh(fac*yh0(iy))*divy
    enddo
  enddo
   !Invert using a full transform in x:
  call revfft(ny+1,nx,xori,xtrig,xfactors)
  call revfft(ny+1,nx,yori,xtrig,xfactors)

   !Add on the x,y coordinate to the deviations:
  do ix=0,nxm1
    do iy=0,ny
      xori(iy,ix)=xori(iy,ix)+glx*dble(ix)-hlx
      yori(iy,ix)=yori(iy,ix)+gly*dble(iy)+ymin
    enddo
  enddo

   !Write coordinates for graphical use:
  open(20,file='coords.r8',form='unformatted', &
      & access='direct',status='replace',recl=2*nbytes)
  write(20,rec=1) 0.d0,xori
  write(20,rec=2) 0.d0,yori
  close(20)

  ! Write conformal domain length for use by the flow_setup script:
  write(*,'(1x,f17.14)') ellx

else

  write(*,*) ' *** Not converging! ***'

endif

contains 

!=======================================================================

double precision function hp(s)
 !Returns the bottom topography height hp at position s:

double precision:: s

hp=a*exp(-b*(s-c)**2)
 !The parameters a, b & c are available from the main calling routine

return
end function

!=======================================================================
      
end program
