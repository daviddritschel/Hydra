program zonal
!  -------------------------------------------------------------------------
!  |   For modes or layers, this routine computes the zonal mean zonal     |
!  |   velocity u, <u>, eddy kinetic energy <(u-<u>)^2+v^2>/2, zonal mean  |
!  |   PV, <q>, eddy enstrophy <(q-<q>)^2>/2, and PV flux <v(q-<q>)> as a  |
!  |   function of y from the PV fields in qq1.r4 & qq2.r4 (for layers 1   |
!  |   and 2 respectively).                                                |
!  |                                                                       |
!  |   Output is to the unformatted, direct-access file "avg.r4" which     |
!  |   contains single-precision numerical data of the form                |
!  |           t, u1, k1, q1, z1, f1, u2, k2, q2, z2, f2                   |
!  |   with each record corresponding to a given time read from the input  |
!  |   files qq1.r4 & qq2.r4, and with u1, etc... being real (single-      |
!  |   precision) arrays dimensioned 0:ny (for all grid lines in y).       |
!  |                                                                       |
!  |   The file ext-zonal.asc max abs values of all of the fields, i.e.    |
!  |   |u1|_max, |k1|_max, ....                                            |
!  -------------------------------------------------------------------------

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use spectral

implicit none
double precision:: fhb(0:ny,0:nxm1)
double precision:: dum
integer:: ix,iy

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
  open(12,file='topo.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) dum,fhb
  close(12)
else
  do ix=0,nxm1
    do iy=0,ny
      fhb(iy,ix)=zero
    enddo
  enddo
endif

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(0:ny,0:nxm1,nz),pp(0:ny,0:nxm1,nz)
double precision:: uu(0:ny,0:nxm1,nz),vv(0:ny,0:nxm1,nz)
double precision:: qq1,qq2,t

 !Diagnostic quantities:
real:: q1r4(0:ny,0:nxm1),q2r4(0:ny,0:nxm1)
real::  u1(0:ny,0:nxm1), v1(0:ny,0:nxm1), q1(0:ny,0:nxm1)
real::  u2(0:ny,0:nxm1), v2(0:ny,0:nxm1), q2(0:ny,0:nxm1)
real:: u1avg(0:ny),k1avg(0:ny),q1avg(0:ny),z1avg(0:ny),f1avg(0:ny)
real:: u2avg(0:ny),k2avg(0:ny),q2avg(0:ny),z2avg(0:ny),f2avg(0:ny)
real:: u1tavg(0:ny),k1tavg(0:ny),q1tavg(0:ny),z1tavg(0:ny),f1tavg(0:ny)
real:: u2tavg(0:ny),k2tavg(0:ny),q2tavg(0:ny),z2tavg(0:ny),f2tavg(0:ny)
real:: yg(0:ny)
real:: u1max,k1max,q1max,z1max,f1max
real:: u2max,k2max,q2max,z2max,f2max
real:: mom1,mom2,momtot
real:: u1bar,u2bar,u1sq,z1sq,u1z1,u2sq,z2sq,u2z2
real:: du1,du2,corr1,corr2
real:: t1,t2,tr4,zfac,hzfac,avgfac

integer:: loop,iread,ix,iy,iopt,ntavg
logical:: timeavg

!---------------------------------------------------------------
 !Open files containing PV field in each layer:
open(31,file='evolution/qq1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(32,file='evolution/qq2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

write(*,*) ' Analyse (1) modes or (2) layers?'
read(*,*) iopt

write(*,*) ' Time period to carry out time average [t1,t2]?'
read(*,*) t1,t2

 !Open output files:
open(22,file='diagnostics/avg.r4',form='unformatted', access='direct', &
                    status='replace',recl=4*(1+10*nyp1))
open(66,file='diagnostics/ext-zonal.asc',status='replace')
open(77,file='diagnostics/momentum.asc',status='replace')
open(88,file='diagnostics/ens-u-correl.asc',status='replace')

timeavg=t2 .ge. t1
if (timeavg) then
   !Initialise time-averaged and zonal-averaged fields:
  do iy=0,ny
    u1tavg(iy)=0.
    k1tavg(iy)=0.
    q1tavg(iy)=0.
    z1tavg(iy)=0.
    f1tavg(iy)=0.
    u2tavg(iy)=0.
    k2tavg(iy)=0.
    q2tavg(iy)=0.
    z2tavg(iy)=0.
    f2tavg(iy)=0.
  enddo
  ntavg=0
  t1=t1-0.001d0*tgsave
  t2=t2+0.001d0*tgsave
endif

 !Frequently used constants:
zfac=1./nx
hzfac=zfac/2.

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,q1r4
  if (iread .ne. 0) exit 
  if (tr4 .gt. t2) exit 
  read(32,rec=loop,iostat=iread) tr4,q2r4
  t=dble(tr4)
  write(*,'(a,f13.5)') ' Processing t = ',t

  do ix=0,nxm1
    do iy=0,ny
      qq(iy,ix,1)=dble(q1r4(iy,ix))
      qq(iy,ix,2)=dble(q2r4(iy,ix))
    enddo
  enddo

   !Invert PV to get velocity field (uu,vv) and streamfunction (pp):
  call main_invert(qq,fhb,t,uu,vv,pp)

  if (iopt .eq. 1) then
     !Compute mode 1 and mode 2 components of u, v & q:
    do ix=0,nxm1
      do iy=0,ny
        u1(iy,ix)=real(vec11*uu(iy,ix,1)+vec12*uu(iy,ix,2))
        u2(iy,ix)=real(vec21*uu(iy,ix,1)+vec22*uu(iy,ix,2))
        v1(iy,ix)=real(vec11*vv(iy,ix,1)+vec12*vv(iy,ix,2))
        v2(iy,ix)=real(vec21*vv(iy,ix,1)+vec22*vv(iy,ix,2))
        q1(iy,ix)=real(vec11*qq(iy,ix,1)+vec12*qq(iy,ix,2))
        q2(iy,ix)=real(vec21*qq(iy,ix,1)+vec22*qq(iy,ix,2))
      enddo
    enddo
  else
     !Just work with layer quantities:
    do ix=0,nxm1
      do iy=0,ny
        u1(iy,ix)=real(uu(iy,ix,1))
        u2(iy,ix)=real(uu(iy,ix,2))
        v1(iy,ix)=real(vv(iy,ix,1))
        v2(iy,ix)=real(vv(iy,ix,2))
        q1(iy,ix)=real(qq(iy,ix,1))
        q2(iy,ix)=real(qq(iy,ix,2))
      enddo
    enddo
  endif

   !Obtain zonal averages; first initialise:
  do iy=0,ny
    u1avg(iy)=0.
    k1avg(iy)=0.
    q1avg(iy)=0.
    z1avg(iy)=0.
    f1avg(iy)=0.
    u2avg(iy)=0.
    k2avg(iy)=0.
    q2avg(iy)=0.
    z2avg(iy)=0.
    f2avg(iy)=0.
  enddo
  
  do ix=0,nxm1
    do iy=0,ny
      u1avg(iy)=u1avg(iy)+u1(iy,ix)
      q1avg(iy)=q1avg(iy)+q1(iy,ix)
      u2avg(iy)=u2avg(iy)+u2(iy,ix)
      q2avg(iy)=q2avg(iy)+q2(iy,ix)
    enddo
  enddo

  do iy=0,ny
    u1avg(iy)=zfac*u1avg(iy)
    q1avg(iy)=zfac*q1avg(iy)
    u2avg(iy)=zfac*u2avg(iy)
    q2avg(iy)=zfac*q2avg(iy)
  enddo

  do ix=0,nxm1
    do iy=0,ny
      k1avg(iy)=k1avg(iy)+(u1(iy,ix)-u1avg(iy))**2+v1(iy,ix)**2
      z1avg(iy)=z1avg(iy)+(q1(iy,ix)-q1avg(iy))**2
      f1avg(iy)=f1avg(iy)+(q1(iy,ix)-q1avg(iy))*v1(iy,ix)
      k2avg(iy)=k2avg(iy)+(u2(iy,ix)-u2avg(iy))**2+v2(iy,ix)**2
      z2avg(iy)=z2avg(iy)+(q2(iy,ix)-q2avg(iy))**2
      f2avg(iy)=f2avg(iy)+(q2(iy,ix)-q2avg(iy))*v2(iy,ix)
    enddo
  enddo

  do iy=0,ny
    k1avg(iy)=hzfac*k1avg(iy)
    z1avg(iy)=hzfac*z1avg(iy)
    f1avg(iy)= zfac*f1avg(iy)
    k2avg(iy)=hzfac*k2avg(iy)
    z2avg(iy)=hzfac*z2avg(iy)
    f2avg(iy)= zfac*f2avg(iy)
  enddo

   !Write diagnostic data for this time:
  write(22,rec=loop) tr4,u1avg,k1avg,q1avg,z1avg,f1avg, &
                         u2avg,k2avg,q2avg,z2avg,f2avg

   !Find max abs values and write to a file:
  u1max=abs(u1avg(0))
  k1max=abs(k1avg(0))
  q1max=abs(q1avg(0))
  z1max=abs(z1avg(0))
  f1max=abs(f1avg(0))
  u2max=abs(u2avg(0))
  k2max=abs(k2avg(0))
  q2max=abs(q2avg(0))
  z2max=abs(z2avg(0))
  f2max=abs(f2avg(0))
  do iy=1,ny
    u1max=max(u1max,abs(u1avg(iy)))
    k1max=max(k1max,abs(k1avg(iy)))
    q1max=max(q1max,abs(q1avg(iy)))
    z1max=max(z1max,abs(z1avg(iy)))
    f1max=max(f1max,abs(f1avg(iy)))
    u2max=max(u2max,abs(u2avg(iy)))
    k2max=max(k2max,abs(k2avg(iy)))
    q2max=max(q2max,abs(q2avg(iy)))
    z2max=max(z2max,abs(z2avg(iy)))
    f2max=max(f2max,abs(f2avg(iy)))
  enddo
  write(66,'(f8.2,10(1x,f11.7))') & 
     tr4,u1max,k1max,q1max,z1max,f1max,u2max,k2max,q2max,z2max,f2max

   !Compute momenta:
  mom1=f12*(u1avg(0)+u1avg(ny))
  mom2=f12*(u2avg(0)+u2avg(ny))
  do iy=1,nym1
    mom1=mom1+u1avg(iy)
    mom2=mom2+u2avg(iy)
  enddo
  mom1=mom1/float(ny)
  mom2=mom2/float(ny)
  if (iopt .eq. 1) then
     !From modal momenta, compute global momentum:
    momtot=h1*(vect11*mom1+vect12*mom2)+alpha*h2*(vect21*mom1+vect22*mom2)
  else
     !From layer momenta, compute global momentum:
    momtot=h1*mom1+alpha*h2*mom2
  endif
  write(77,'(f8.2,3(1x,f11.7))') tr4,mom1,mom2,momtot

   !Correlate the zonal velocity (minus its mean over y) with vorticity:
  u1bar=f12*(u1avg(0)+u1avg(ny))
  u2bar=f12*(u2avg(0)+u2avg(ny))
  do iy=1,ny-1
    u1bar=u1bar+u1avg(iy)
    u2bar=u2bar+u2avg(iy)
  enddo
  u1bar=u1bar/dble(ny)
  u2bar=u2bar/dble(ny)
  u1sq=zero
  z1sq=zero
  u1z1=zero
  u2sq=zero
  z2sq=zero
  u2z2=zero
  do iy=1,ny-1
    du1=u1avg(iy)-u1bar
    u1sq=u1sq+du1**2
    z1sq=z1sq+z1avg(iy)**2
    u1z1=u1z1+du1*z1avg(iy)
    du2=u2avg(iy)-u2bar
    u2sq=u2sq+du2**2
    z2sq=z2sq+z2avg(iy)**2
    u2z2=u2z2+du2*z2avg(iy)
  enddo
  corr1=u1z1/sqrt(u1sq*z1sq)
  corr2=u2z2/sqrt(u2sq*z2sq)
  write(88,'(f8.2,2(1x,f11.7))') tr4,corr1,corr2

  if (timeavg) then
     !Accumulate time average if t1-eps < t < t2+eps where eps = tgsave/1000:
    if (t1 .lt. t .and. t .lt. t2) then
      ntavg=ntavg+1
      do iy=0,ny
        u1tavg(iy)=u1tavg(iy)+u1avg(iy)
        k1tavg(iy)=k1tavg(iy)+k1avg(iy)
        q1tavg(iy)=q1tavg(iy)+q1avg(iy)
        z1tavg(iy)=z1tavg(iy)+z1avg(iy)
        f1tavg(iy)=f1tavg(iy)+f1avg(iy)
        u2tavg(iy)=u2tavg(iy)+u2avg(iy)
        k2tavg(iy)=k2tavg(iy)+k2avg(iy)
        q2tavg(iy)=q2tavg(iy)+q2avg(iy)
        z2tavg(iy)=z2tavg(iy)+z2avg(iy)
        f2tavg(iy)=f2tavg(iy)+f2avg(iy)
      enddo
    endif
  endif
enddo
 !End of loop over time

if (timeavg) then
   !Finalise time-averaged and zonal-averaged fields:
  avgfac=one/dble(ntavg)
  do iy=0,ny
    yg(iy)=gly*dble(iy)/elly-f12
    u1tavg(iy)=u1tavg(iy)*avgfac
    k1tavg(iy)=k1tavg(iy)*avgfac
    q1tavg(iy)=q1tavg(iy)*avgfac
    z1tavg(iy)=z1tavg(iy)*avgfac
    f1tavg(iy)=f1tavg(iy)*avgfac
    u2tavg(iy)=u2tavg(iy)*avgfac
    k2tavg(iy)=k2tavg(iy)*avgfac
    q2tavg(iy)=q2tavg(iy)*avgfac
    z2tavg(iy)=z2tavg(iy)*avgfac
    f2tavg(iy)=f2tavg(iy)*avgfac
  enddo

   !Write data for use by plotcol, with y in the second row:
  open(41,file='diagnostics/u1avg.asc',status='replace')
  open(42,file='diagnostics/k1avg.asc',status='replace')
  open(43,file='diagnostics/q1avg.asc',status='replace')
  open(44,file='diagnostics/z1avg.asc',status='replace')
  open(45,file='diagnostics/f1avg.asc',status='replace')
  open(51,file='diagnostics/u2avg.asc',status='replace')
  open(52,file='diagnostics/k2avg.asc',status='replace')
  open(53,file='diagnostics/q2avg.asc',status='replace')
  open(54,file='diagnostics/z2avg.asc',status='replace')
  open(55,file='diagnostics/f2avg.asc',status='replace')
  do iy=0,ny
    write(41,'(f14.9,1x,f14.9)') u1tavg(iy),yg(iy)
    write(42,'(f14.9,1x,f14.9)') k1tavg(iy),yg(iy)
    write(43,'(f14.9,1x,f14.9)') q1tavg(iy),yg(iy)
    write(44,'(f14.9,1x,f14.9)') z1tavg(iy),yg(iy)
    write(45,'(f14.9,1x,f14.9)') f1tavg(iy),yg(iy)
    write(51,'(f14.9,1x,f14.9)') u2tavg(iy),yg(iy)
    write(52,'(f14.9,1x,f14.9)') k2tavg(iy),yg(iy)
    write(53,'(f14.9,1x,f14.9)') q2tavg(iy),yg(iy)
    write(54,'(f14.9,1x,f14.9)') z2tavg(iy),yg(iy)
    write(55,'(f14.9,1x,f14.9)') f2tavg(iy),yg(iy)
  enddo
  close(41)
  close(42)
  close(43)
  close(44)
  close(45)
  close(51)
  close(52)
  close(53)
  close(54)
  close(55)
endif

 !Close output files:
close(22)
close(66)
close(77)
close(88)

write(*,*)
write(*,*) ' The results are ready in avg.r4; to view the results type'
write(*,*) ' zv'
write(*,*)
write(*,*) ' The extremal values vs time are listed in ext-zonal.asc'

if (timeavg) then
  write(*,*)
  write(*,*) ' Time averaged profiles are in u1avg.asc, k1avg.asc, etc.'
endif

write(*,*)
write(*,*) ' The layer/modal momenta and global momentum is in momentum.asc'
write(*,*) ' and the enstrophy-zonal flow correlation is in ens-u-correl.asc'

return
end subroutine

 !End main program
end program
!=======================================================================
