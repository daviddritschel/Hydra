program v3
!-----------------------------------------------------------------
!    Generates three Gaussian vortices of relative vorticity
!    with zero divergence in geostrophic balance.
!-----------------------------------------------------------------

use spectral

implicit double precision(a-h,o-z)
double precision:: qs(ng,ng),qq(ng,ng),zz(ng,ng),hh(ng,ng)

 !Initialise inversion constants and arrays:
call init_spectral

write(*,*) ' Vorticity amplitude/f for vortex 1? '
read(*,*) a1
write(*,*) ' x_1 relative to domain width/18? '
read(*,*) x1
write(*,*) ' y_1 relative to domain width/18? '
read(*,*) y1
write(*,*) ' R_1 relative to domain width/18? '
read(*,*) r1

write(*,*)
write(*,*) ' Vorticity amplitude/f for vortex 2? '
read(*,*) a2
write(*,*) ' x_2 relative to domain width/18? '
read(*,*) x2
write(*,*) ' y_2 relative to domain width/18? '
read(*,*) y2
write(*,*) ' R_2 relative to domain width/18? '
read(*,*) r2

write(*,*)
write(*,*) ' Vorticity amplitude/f for vortex 3? '
read(*,*) a3
write(*,*) ' x_3 relative to domain width/18? '
read(*,*) x3
write(*,*) ' y_3 relative to domain width/18? '
read(*,*) y3
write(*,*) ' R_3 relative to domain width/18? '
read(*,*) r3

! Redefine amplitudes including f:
a1=a1*cof
a2=a2*cof
a3=a3*cof

! Redefine positions and radii including domain width/18:
fac=twopi/18.d0
x1=x1*fac-pi
y1=y1*fac-pi
s1=f12/(r1*fac)**2
x2=x2*fac-pi
y2=y2*fac-pi
s2=f12/(r2*fac)**2
x3=x3*fac-pi
y3=y3*fac-pi
s3=f12/(r3*fac)**2

! Generate vorticity field and enforce periodicity (use zz):
zz=zero
do j=-1,1
  xoff=twopi*dble(j)-pi
  do i=-1,1
    yoff=twopi*dble(i)-pi
    do ix=1,ng
      xg=gl*dble(ix-1)+xoff
      do iy=1,ng
        yg=gl*dble(iy-1)+yoff
        zz(iy,ix)=zz(iy,ix)+a1*exp(-s1*((xg-x1)**2+(yg-y1)**2))+ &
                            a2*exp(-s2*((xg-x2)**2+(yg-y2)**2))+ &
                            a3*exp(-s3*((xg-x3)**2+(yg-y3)**2))
      enddo
    enddo
  enddo
enddo

! Remove mean value:
qbar=dsumi*sum(zz)
zz=zz-qbar

! Copy to qq then transform to spectral space as qs:
qq=zz
call ptospc(ng,ng,qq,qs,xfactors,yfactors,xtrig,ytrig)

! Obtain height anomaly by geostrophic balance:
geo=cof/cgw**2
qs=geo*rlap*qs

! Return to physical space as hh:
call spctop(ng,ng,qs,hh,xfactors,yfactors,xtrig,ytrig)

! Define the PV anomaly field:
qq=(zz+cof)/(one+hh)-cof

! Write PV:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,qq
close(11)

! Write zero divergence and acceleration divergence (geostrophic balance):
zz=zero
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,zz
close(11)

open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,zz
close(11)

end program
