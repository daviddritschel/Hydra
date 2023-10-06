module evolution

! Module contains subroutines to evolve PV according to the 
! algorithm detailed in casl.f90.

use common

implicit none

 !Energies:
double precision:: ekepre,ekepost,apepre,apepost

 !Physical fields:
double precision::    qq(0:ny,0:nxm1),qqpre(0:ny,0:nxm1)
double precision:: qspre(0:ny,0:nxm1),qdpre(0:ny,0:nxm1)
double precision::    uu(0:ny,0:nxm1),   vv(0:ny,0:nxm1)
double precision::    zz(0:ny,0:nxm1),zzpre(0:ny,0:nxm1)
double precision::    pp(0:ny,0:nxm1),   qc(0:ny,0:nxm1)

 !Semi-Lagrangian advection arrays:
double precision::  x0(0:ny,0:nxm1), y0(0:ny,0:nxm1)
double precision:: ula(0:ny,0:nxm1),vla(0:ny,0:nxm1)

!Internal subroutine definitions (inherit global variables):

contains 

!=============================================================
subroutine advect

! Main subroutine for advecting fields and contours

implicit none

 !Local variables:
double precision,parameter:: twistmax=2.5d0
!      twistmax: the maximum value of the time integral of |zeta|_max
!                between regularisations of the contours.
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops 
!      to rebuild the contours in a separate memory space.
integer:: ireg,igsave,icsave

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Used for regularising contours:
twist=zero

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Perform surgery & field reset 
  if (twist .gt. twistmax) then
    ireg=ireg+1
     !Don't continue if maximum number of regularisations reached;
     !it is time to recontour (pass control back to main program):
    if (ireg .eq. nregmax) then
       !Prepare for recontouring:
      qd=qq-qc
      qs=qq
       !Exit module and go to recontouring:
      return
    endif

     !Regularise contours (surgery + node redistribution):
    call surgery
     !Record contour complexity to complexity.asc:
    write(14,'(1x,f13.5,1x,i8,1x,i9)') t,nq,nptq

     !Convert PV contours to gridded values (qc):
    call con2grid(qc)
     !Reset qs and qd:
    call reset(qc,qs,qd,qavg)

     !Copy gridded PV fields to old time level:
    call copyfields

     !Update twist parameter:
    twist=twist-twistmax
  endif

   !Adjust timestep (dt) on maximum vorticity magnitude:
  call adapt(igsave,icsave)
   !Advect PV from time t to t + dt:
  call advance

   !Update the time:
  t=t+dt

   !Possibly save PV & energy at chosen save time (tgrid):
  if (igsave .eq. 1) call savegrid

   !Possibly save contours and residual PV (qd) for post processing:
  if (icsave .eq. 1) call savecont
  
   !Copy new fields into previous time:
  call copyfields
enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
return
end subroutine

!=======================================================================

subroutine copyfields

! Copies new fields into previous time

implicit none

qqpre=qq
qdpre=qd
qspre=qs
zzpre=zz

return
end subroutine

!=======================================================================

subroutine init

! Initialises quantities needed for normal time integration following 
! contour regeneration

implicit none

!---------------------------------------------
 !Record contour complexity to complexity.asc:
write(14,'(1x,f13.5,1x,i8,1x,i9)') t,nq,nptq

!---------------------------------------------------
 !Convert PV contours to gridded values (qc):
call con2grid(qc)

 !Define residual PV qd = qs-qc-F[qs-qc]
ula=qs-qc
vla=ula

call filter(ula,0,2)

 !Finish qd calculation & copy fields to old time level:
qd=vla-ula
qqpre=qs
qdpre=qd
qspre=qs

!--------------------------------------------------------
 !Get the initial velocity (uu,vv) and streamfunction pp:
call inversion
zzpre=zz

return
end subroutine

!=======================================================================

subroutine inversion

! Combines the contour and grid PV fields then inverts the PV to
! obtain the velocity and streamfunction.

implicit none

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qc)

 !Combine fields to update qq with full field:
call combine(qq,qc,qs,qd,qavg)

 !Invert PV to obtain velocity field (uu,vv), streamfunction (pp)
 !and relative vorticity (zz):
call main_invert(qq,uu,vv,pp,zz)

return
end subroutine

!=======================================================================
      
subroutine advance

! Advances PV contours and gridded fields to time t+dt by contour 
! advection and by trajectory integration using a standard 
! semi-Lagrangian (SL) scheme.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
double precision:: uq(nptq),vq(nptq),xqm(nptq),yqm(nptq)

!------------------------------------------------------------------------
 !Increments in grid units needed below for trajectory integration:
gcx=dt*glxi
gcy=dt*glyi
hgcx=f12*gcx
hgcy=f12*gcy

 !Copy current velocity field into (ula,vla) for use below then compute
 !Euler backward predictor for departure grid point (x0,y0):
do ix=0,nxm1
  do iy=0,ny
    ula(iy,ix)=uu(iy,ix)
    vla(iy,ix)=vv(iy,ix)
    x0(iy,ix)=mod(xigmax + xig(ix)-gcx*uu(iy,ix),xigmax)
    y0(iy,ix)=max(zero,min(yig(iy)-gcy*vv(iy,ix),yigmax))
  enddo
enddo
 !Note, (uu,vv) is used since we have no other velocity field available

 !Prepare contour evolution; get velocity on PV contour nodes:
if (nptq .gt. 0) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xq(i)+hfdt*uq(i)
    xqm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqm(i)=min(ymax,max(ymin,yq(i)+hfdt*vq(i)))
    xx=xq(i)+dt*uq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=min(ymax,max(ymin,yq(i)+dt*vq(i)))
  enddo
endif

 !Iterate to converge on implicit trapezoidal integration:
do iter=1,niter
   !Obtain qs & qd at time t + dt:
  call sl_step

   !Obtain qq & hence uu & vv at time t + dt:
  call inversion

   !Correct departure grid point (x0,y0):
  do ix=0,nxm1
    do iy=0,ny
       !Obtain old velocity (time t) at the departure grid point using
       !bi-linear interpolation of (ula,vla):
      ix0=int(x0(iy,ix))
      ix1=ixp(ix0)
      px=x0(iy,ix)-dble(ix0)
      pxc=one-px

      iy0=int(y0(iy,ix))
      iy1=iyp(iy0)
      py=y0(iy,ix)-dble(iy0)
      pyc=one-py

      uod=pyc*(pxc*ula(iy0,ix0)+px*ula(iy0,ix1)) &
          +py*(pxc*ula(iy1,ix0)+px*ula(iy1,ix1))

      vod=pyc*(pxc*vla(iy0,ix0)+px*vla(iy0,ix1)) &
          +py*(pxc*vla(iy1,ix0)+px*vla(iy1,ix1))

      x0(iy,ix)=mod(xigmax + xig(ix)-hgcx*(uod+uu(iy,ix)),xigmax)
      y0(iy,ix)=max(zero,min(yig(iy)-hgcy*(vod+vv(iy,ix)),yigmax))
       !(uu,vv) is the new velocity (time t+dt) at the arrival grid point
    enddo
  enddo

   !Update the PV contour points:
  if (nptq .gt. 0) then
    call velint(uu,vv,uq,vq)
    do i=1,nptq
      xx=xqm(i)+hfdt*uq(i)
      xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
      yq(i)=min(ymax,max(ymin,yqm(i)+hfdt*vq(i)))
    enddo
  endif

enddo

 !Obtain final corrected qs & qd at time t + dt:
call sl_step

 !Obtain final corrected uu & vv at time t + dt from qq:
call inversion

 !Update the PV contour points:
if (nptq .gt. 0) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqm(i)+hfdt*uq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=min(ymax,max(ymin,yqm(i)+hfdt*vq(i)))
  enddo
endif

 !Call con2grid to get updated contour PV (qc):
call con2grid(qc)

 !Combine fields to update qq with full field:
call combine(qq,qc,qs,qd,qavg)

return
end subroutine

!=======================================================================
      
subroutine sl_step

! Interpolates qs & qd at points (x0,y0) to perform SL integration
! from t to t+dt for qs & qd.

implicit none

!-------------------------------------------------------------------------
 !Integrate qs & qd using bi-cubic Lagrange interpolation of qspre & qdpre
 !at x0,y0:
call interpol(qspre,qs,x0,y0)
call interpol(qdpre,qd,x0,y0)
 !Here, we obtain only the adiabatic evolution

return
end subroutine

!=======================================================================

subroutine adapt(igsave,icsave)

! Adapts the time step to ensure dt < dtfac/max(|zeta|_max)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local parameter used for setting time step:
double precision,parameter:: dtfac=pi/40.d0
integer:: ix,iy

!------------------------------------------------------------------------
 !Compute max abs and rms vorticity:
zzmax=small
zzrms=zero
do ix=0,nxm1
  do iy=0,nym1
    zzmax=max(zzmax,abs(zz(iy,ix)))
    zzrms=zzrms+zz(iy,ix)**2
  enddo
enddo
zzrms=sqrt(dsumi*zzrms)

 !Time step needed for accuracy:
dtacc=dtfac/max(zzmax,srwfm)
 !The restriction on the maximum Rossby wave frequency (srwfm)
 !ensures that the fastest Rossby wave frequency is resolved.

 !Choose a new time step, limiting it to dtmax in parameters.f90:
dt=min(dtacc,dtmax)
hfdt=dt/two

 !Increment the integral of max|zz|:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
 !Record various diagnostics to monitor.asc:
write(12,'(1x,f13.5,2(1x,f14.8))') t,zzmax,zzrms

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
tgrid=tinit+tgsave*(dble(igrec)+1.d-10)
 !A small number is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  igsave=1
   !Compute kinetic & potential energy:
  call energy(ekepre,apepre)
else
   !Do not save data:
  igsave=0
endif

 !Set flag to save contour data every tcsave time units:
tcont=tinit+tcsave*dble(icrec)
if (t .le. tcont .and. t+dt .gt. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  icsave=1
else
   !Do not save data:
  icsave=0
endif

return
end subroutine

!=======================================================================

subroutine energy(eke,ape)

! Computes the kinetic and potential energy from the velocity field
! (uu,vv) and the streamfunction pp.

! Note: the energies are normalised by rho*H.

implicit none

double precision:: eke,ape,ul2,vl2,pl2
integer:: ix,iy

 !Kinetic energy:
call l2norm(uu,ul2)
call l2norm(vv,vl2)
eke=f12*(ul2+vl2)

 !Potential energy:
call l2norm(pp,pl2)
ape=f12*kdsq*pl2

return
end subroutine

!=======================================================================

subroutine savegrid

! Saves various quantities at the desired save time to files 

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: qspec(0:max(nx,ny))
real:: qr4(0:ny,0:nxm1),tr4

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Needed for real*4 writes below:
tr4=real(tgrid)

 !Record to write:
irec=igrec+1

!-----------------------------------------------------
 !Store PV:
qr4=real(pt*qqpre+ptc*qq)
write(31,rec=irec) tr4,qr4

 !Store relative vorticity:
qr4=real(pt*zzpre+ptc*zz)
write(32,rec=irec) tr4,qr4

!----------------------------------------------------------
 !Compute kinetic & potential energy:
call energy(ekepost,apepost)

 !Compute time interpolated energies:
apot=pt*apepre+ptc*apepost
ekin=pt*ekepre+ptc*ekepost

 !Write diagnostics to ene.asc:
write(15,'(f8.2,3(1x,f14.9))') tgrid,ekin,apot,ekin+apot

write(*,'(a,f8.2,3(a,f11.7))') &
    & ' t = ',tgrid,'  K = ',ekin,'  P = ',apot,'  E = K+P = ',ekin+apot

!-----------------------------------------------------
 !Compute PV anomaly and its spectrum:
do ix=0,nxm1
  do iy=0,ny
    vla(iy,ix)=qq(iy,ix)-bety(iy)
  enddo
enddo

 !Compute the 1d PV spectrum:
call spec1d_fc(vla,qspec)

sqspec=zero
do k=0,kmax
  sqspec=sqspec+qspec(k)
   !Normalise to take into account uneven sampling of wavenumbers
   !in each shell [k-1/2,k+1/2]:
  qspec(k)=spmf(k)*qspec(k)
enddo
sqspec=8.d0*sqspec*dsumi
 !Write out spectrum to file:
write(51,'(f8.2,1x,e14.7,1x,i5)') tgrid,sqspec,kmaxred 
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),log10(qspec(k))
enddo
 !Note: alk(k) = log_10(k)

!--------------------------------------------------------------
 !Increment counter for direct file access:
igrec=igrec+1

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves PV contours and residual PV for use by congen.f90 & diagnostics

implicit double precision(a-h,o-z)
implicit integer(i-n)

real:: qdr4(0:ny,0:nxm1),tr4
integer:: iopq(nq)
character(len=3):: pind

write(*,'(a,f13.5)') ' Saving contours at t = ',t

write(pind(1:3),'(i3.3)') icrec

tr4=real(t)

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f13.5,2(1x,f16.12))') nq,nptq,t,dq,qavg

 !Save PV contours if any exist:
if (nq .gt. 0) then
   !First form iopq; open/closed indicator:
  do j=1,nq
    iopq(j)=nextq(i2q(j))/i1q(j)
     !iopq = 0 for an open contour, and 1 for a closed one
  enddo
  open(81,file='contours/qqindex'//pind,form='unformatted', &
        access='direct',status='replace',recl=20*nq)
  write(81,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq),iopq(1:nq)
  close(81)

  open(82,file='contours/qqnodes'//pind,form='unformatted', &
        access='direct',status='replace',recl=16*nptq)
  write(82,rec=1) xq(1:nptq),yq(1:nptq)
  close(82)
endif

 !Save residual needed to build ultra-fine-grid PV with congen:
open(83,file='contours/qqresi'//pind,form='unformatted', &
      access='direct',status='replace',recl=nbytes)
qdr4=real(qq-qc)
write(83,rec=1) tr4,qdr4
close(83)

!--------------------------------------------------------------
 !Increment counter for naming next direct-access output files:
icrec=icrec+1

return
end subroutine

!=======================================================================

 !Main end module
end module
