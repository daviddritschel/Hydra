module evolution

! Module contains subroutines to evolve bb and zz fields according to the 
! algorithm detailed in casl.f90.

use common

implicit none

 !Energies:
double precision:: ekepre,ekepost,apepre,apepost,dkdtpre,dkdtpost

 !Physical fields:
double precision:: zz(0:ny,0:nxm1),zc(0:ny,0:nxm1),sz(0:ny,0:nxm1)
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: zzpre(0:ny,0:nxm1),zspre(0:ny,0:nxm1),zdpre(0:ny,0:nxm1)
double precision:: bbpre(0:ny,0:nxm1),szpre(0:ny,0:nxm1)

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
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops 
!      to rebuild the contours in a separate memory space.
integer:: istep,ireg,igsave,icsave,ix,iy
double precision:: zzl1,zzl2

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Counter used for resetting fields and performing contour surgery:
istep=0

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

   !Perform surgery & field reset every nstep time steps:
  istep=mod(istep,nstep)+1
  if (istep .eq. nstep) then
    ireg=ireg+1
     !Don't continue if maximum number of regularisations reached;
     !it is time to recontour (pass control back to main program):
    if (ireg .eq. nregmax) then
       !Update average vorticity:
      call average(zz,zavg)
       !Prepare vorticity residual for recontouring:
      do ix=0,nxm1
        do iy=0,ny
          zd(iy,ix)=zz(iy,ix)-zc(iy,ix)
          zs(iy,ix)=zz(iy,ix)
        enddo
      enddo
       !Compute contour interval for vorticity:
      call l1norm(zz,zzl1)
      call l2norm(zz,zzl2)
      zjump=(zzl2/zzl1)/dble(ncontz)
       !ncontz is set in parameters.f90. 
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the buoyancy contours (surgery + node redistribution):
    call surgery(xb,yb,nextb,indb,npb,i1b,i2b,nb,nptb)
    call surgery(xz,yz,nextz,indz,npz,i1z,i2z,nz,nptz)
     !Record contour complexity to complexity.dat:
    write(14,'(1x,f12.5,2(1x,i8,1x,i9))') t,nb,nptb,nz,nptz

     !Convert vorticity contours to gridded values (zc):
    call con2grid(zc,xz,yz,zjump,zavg,nextz,nptz,-1)
     !Reset zs and zd:
    call reset(zc,zs,zd,zavg)

     !Copy gridded vorticity fields to old time level:
    do ix=0,nxm1
      do iy=0,ny
        zzpre(iy,ix)=zs(iy,ix) 
        zdpre(iy,ix)=zd(iy,ix) 
        zspre(iy,ix)=zs(iy,ix) 
      enddo
    enddo
  endif

   !Adjust timestep (dt) on maximum vorticity magnitude:
  call adapt(igsave,icsave)
   !Advect buoyancy & vorticity from time t to t + dt:
  call advance

   !Update the time:
  t=t+dt

   !Possibly save buoyancy, vorticity & energy at chosen save time (tgrid):
  if (igsave .eq. 1) call savegrid

   !Possibly save contours and residual vorticity (zd) for post processing:
  if (icsave .eq. 1) call savecont
  
   !Copy new fields into previous time:
  do ix=0,nxm1
    do iy=0,ny
      szpre(iy,ix)=sz(iy,ix)
      zzpre(iy,ix)=zz(iy,ix)
      zdpre(iy,ix)=zd(iy,ix)
      zspre(iy,ix)=zs(iy,ix)
    enddo
  enddo
enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
return
end subroutine

!=======================================================================

subroutine init

! Initialises quantities needed for normal time integration following 
! contour regeneration

implicit double precision(a-h,o-z)
implicit integer(i-n)

!---------------------------------------------
 !Record contour complexity to complexity.dat:
write(14,'(1x,f12.5,2(1x,i8,1x,i9))') t,nb,nptb,nz,nptz

!---------------------------------------------------
 !Convert vorticity contours to gridded values (zc):
call con2grid(zc,xz,yz,zjump,zavg,nextz,nptz,-1)

 !Define residual vorticity zd = zs-zc-F[zs-zc]
do ix=0,nxm1
  do iy=0,ny
    ula(iy,ix)=zs(iy,ix)-zc(iy,ix)
    vla(iy,ix)=ula(iy,ix)
  enddo
enddo

call filter(ula,0,2)

do ix=0,nxm1
  do iy=0,ny
    zd(iy,ix)=vla(iy,ix)-ula(iy,ix)
  enddo
enddo

!-----------------------------------------------------
 !Copy gridded vorticity to old time level:
do ix=0,nxm1
  do iy=0,ny
    zzpre(iy,ix)=zs(iy,ix) 
    zdpre(iy,ix)=zd(iy,ix) 
    zspre(iy,ix)=zs(iy,ix) 
  enddo
enddo

!------------------------------------------------------
 !Calculate the source term for zeta (szpre):
call getzzsrc(szpre,t)
 !Note: zzpre and szpre are needed by subroutine advance.

 !Get the initial velocity field (uu,vv):
if (zjump .gt. zero) then 
  call inversion
else
   !Here vorticity is zero identically - no need to do inversion:
  do ix=0,nxm1
    do iy=0,ny
      uu(iy,ix)=zero
      vv(iy,ix)=zero 
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine inversion

! Inverts Laplace's operator on vorticity (zz) to obtain the 
! streamfunction (pp) and the velocity (uu,vv) = (-dpp/dy,dpp/dx).

implicit double precision(a-h,o-z)
implicit integer(i-n)

!------------------------------------------------------------
 !Call con2grid to get updated contour vorticity (zc):
call con2grid(zc,xz,yz,zjump,zavg,nextz,nptz,-1)

 !Combine fields to update zz with full field:
call combine(zz,zc,zs,zd,zavg)

 !Invert Laplace operator on zz spectrally and obtain velocity field: 
call main_invert(zz,uu,vv)      

return
end subroutine

!=======================================================================
      
subroutine advance

! Computes bb(t+dt) by contour advection and zz(t+dt) by trajectory 
! integration using a standard semi-Lagrangian (SL) scheme.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
double precision:: ub(nptb),vb(nptb),xbm(nptb),ybm(nptb)
double precision:: uz(nptz),vz(nptz),xzm(nptz),yzm(nptz)

!------------------------------------------------------------------------
 !Copy current velocity field into (ula,vla) for use in subroutines below:
do ix=0,nxm1
  do iy=0,ny
    ula(iy,ix)=uu(iy,ix)
    vla(iy,ix)=vv(iy,ix)
  enddo
enddo

 !Increments in grid units needed below for trajectory integration:
gcx=dt*glxi
gcy=dt*glyi
hgcx=f12*gcx
hgcy=f12*gcy

 !Compute Euler backward predictor for departure grid point (x0,y0):
do ix=0,nxm1
  do iy=0,ny
    x0(iy,ix)=mod(xigmax+xig(ix)-gcx*uu(iy,ix),xigmax)
    y0(iy,ix)=max(zero,min(yig(iy)-gcy*vv(iy,ix),yigmax))
  enddo
enddo
 !Note, (uu,vv) is used since we have no other velocity field available

 !Prepare contour evolution; get velocity on buoyancy contour nodes:
call velint(uu,vv,xb,yb,ub,vb,nptb)
do i=1,nptb
  xx=xb(i)+hfdt*ub(i)
  xbm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  ybm(i)=min(ymax,max(ymin,yb(i)+hfdt*vb(i)))
  xx=xb(i)+dt*ub(i)
  xb(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yb(i)=min(ymax,max(ymin,yb(i)+dt*vb(i)))
enddo

 !Prepare contour evolution; get velocity on vorticity contour nodes:
if (nptz .gt. 0) then
  call velint(uu,vv,xz,yz,uz,vz,nptz)
  do i=1,nptz
    xx=xz(i)+hfdt*uz(i)
    xzm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yzm(i)=min(ymax,max(ymin,yz(i)+hfdt*vz(i)))
    xx=xz(i)+dt*uz(i)
    xz(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yz(i)=min(ymax,max(ymin,yz(i)+dt*vz(i)))
  enddo
endif

 !Iterate to converge on implicit trapezoidal integration:
do iter=1,niter
   !Obtain zs & zd at time t + dt:
  call sl_step

   !Obtain zz & hence uu & vv at time t + dt:
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
      &   +py*(pxc*ula(iy1,ix0)+px*ula(iy1,ix1))

      vod=pyc*(pxc*vla(iy0,ix0)+px*vla(iy0,ix1)) &
      &   +py*(pxc*vla(iy1,ix0)+px*vla(iy1,ix1))

      x0(iy,ix)=mod(xigmax+xig(ix)-hgcx*(uod+uu(iy,ix)),xigmax)
      y0(iy,ix)=max(zero,min(yig(iy)-hgcy*(vod+vv(iy,ix)),yigmax))
       !(uu,vv) is the new velocity (time t+dt) at the arrival grid point
    enddo
  enddo

   !Update the buoyancy contour points:
  call velint(uu,vv,xb,yb,ub,vb,nptb)
  do i=1,nptb
    xx=xbm(i)+hfdt*ub(i)
    xb(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yb(i)=min(ymax,max(ymin,ybm(i)+hfdt*vb(i)))
  enddo

   !Update the vorticity contour points:
  if (nptz .gt. 0) then
  call velint(uu,vv,xz,yz,uz,vz,nptz)
  do i=1,nptz
    xx=xzm(i)+hfdt*uz(i)
    xz(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yz(i)=min(ymax,max(ymin,yzm(i)+hfdt*vz(i)))
  enddo
  endif

enddo

 !Obtain final corrected zs & zd at time t + dt:
call sl_step

 !Obtain final corrected uu & vv at time t + dt from zz:
call inversion

 !Update the buoyancy contour points:
call velint(uu,vv,xb,yb,ub,vb,nptb)
do i=1,nptb
  xx=xbm(i)+hfdt*ub(i)
  xb(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yb(i)=min(ymax,max(ymin,ybm(i)+hfdt*vb(i)))
enddo

 !Update the vorticity contour points:
if (nptz .gt. 0) then
  call velint(uu,vv,xz,yz,uz,vz,nptz)
  do i=1,nptz
    xx=xzm(i)+hfdt*uz(i)
    xz(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yz(i)=min(ymax,max(ymin,yzm(i)+hfdt*vz(i)))
  enddo
endif

 !Call con2grid to get updated contour vorticity (zc):
call con2grid(zc,xz,yz,zjump,zavg,nextz,nptz,-1)

 !Combine fields to update zz with full field:
call combine(zz,zc,zs,zd,zavg)

return
end subroutine

!=======================================================================
      
subroutine sl_step

! Interpolates zz at points (x0,y0) and adds the source integral
! from t to t+dt to zz

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: tmp(0:ny,0:nxm1)

 !Integrate zs & zd using bi-cubic Lagrange interpolation of zspre & zdpre
 !at x0,y0:
call interpol(zspre,zs,x0,y0)
call interpol(zdpre,zd,x0,y0)
 !Here, we obtain only the adiabatic evolution; add next the source:

call getzzsrc(sz,t+dt)

do ix=0,nxm1
  do iy=0,ny
     !Obtain old source (at time t) at the departure grid point using
     !bi-linear interpolation of szpre:
    ix0=int(x0(iy,ix))
    ix1=ixp(ix0)
    px=x0(iy,ix)-dble(ix0)
    pxc=one-px

    iy0=int(y0(iy,ix))
    iy1=iyp(iy0)
    py=y0(iy,ix)-dble(iy0)
    pyc=one-py

    szod=pyc*(pxc*szpre(iy0,ix0)+px*szpre(iy0,ix1)) &
      &  +py*(pxc*szpre(iy1,ix0)+px*szpre(iy1,ix1))

     !Integrate in time using trapezoidal rule (sz is the new source
     !at time t+dt) and add to zd:
    zd(iy,ix)=zd(iy,ix)+hfdt*(szod+sz(iy,ix))
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine adapt(igsave,icsave)

! Adapts the time step to ensure dt < dtfac/max(|zeta|_max)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Compute accurate advection time step (dtacc):
zzmax=small
do ix=0,nxm1
  do iy=0,ny
    zzmax=max(zzmax,abs(zz(iy,ix)))
  enddo
enddo
 !Note: szpre contains the vorticity source at the current time
dtacc=dtfac/zzmax

!---------------------------------------------------------------------
 !Choose a new time step: 
if (dt .gt. zero) then 
  dt=min(dtacc,dtmax)
  if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
else
   !Calculate a maximum timestep (dtmax) based on grad(bb):
  call bxderiv(bb,ula)
  call byderiv(bb,vla)
  dbmax=small
  do ix=0,nxm1
    do iy=0,ny
      dbmax=max(dbmax,ula(iy,ix)**2+vla(iy,ix)**2)
    enddo
  enddo 
   !Limit max timestep to a data save time
  dtmax=min(tgsave,tcsave,dtfac/sqrt(sqrt(dbmax)))
  dt=min(dtacc,dtmax)
endif 
hfdt=dt/two

!---------------------------------------------------------------------
 !Record max |z| to monitor.dat:
write(12,'(1x,f12.5,1x,1p,e14.7)') t,zzmax

!---------------------------------------------------------------------

 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
tgrid=tgsave*dble(itime)+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  igsave=1
   !Compute current gridded buoyancy and store in bbpre:
  call con2grid(bbpre,xb,yb,bjump,bavg,nextb,nptb,iene)
   !With iene = 0, this routine also computes the reference state 
   !potential energy; with iene >= 0, the available potential energy
   !is computed and stored in the common variable ape.
  apepre=ape
  iene=1
   !Compute kinetic energy:
  call l2norm(uu,uul2)
  call l2norm(vv,vvl2)
  ekepre=f12*(uul2+vvl2)
   !Compute rate of change of kinetic energy:
  call binorm(vv,bbpre,dkdtpre)
else
   !Do not save data:
  igsave=0
endif

 !Set flag to save contour data every tgsave time units:
itime=int((t+dt)/tcsave)
tcont=tcsave*dble(itime)-small
if (t .lt. tcont .and. t+dt .ge. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  icsave=1
else
   !Do not save data:
  icsave=0
endif

return
end subroutine

!=======================================================================
      
subroutine savegrid

! Saves zz & bb at the desired save time to files (process with image):

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: zspec(0:max(nx,ny))
real:: zzr4(0:ny,0:nxm1),bbr4(0:ny,0:nxm1),tr4

 !Increment counter for direct file access:                                                        
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Compute available potential energy through a call to con2grid:
call con2grid(bb,xb,yb,bjump,bavg,nextb,nptb,1)
apepost=ape
 !Compute kinetic energy
call l2norm(uu,uul2)
call l2norm(vv,vvl2)
ekepost=f12*(uul2+vvl2)

 !Compute time interpolated energies:
apot=pt*apepre+ptc*apepost
ekin=pt*ekepre+ptc*ekepost

 !Compute rate of change of kinetic energy:
call binorm(vv,bb,dkdtpost)
dkdt=pt*dkdtpre+ptc*dkdtpost

 !Store vorticity and buoyancy at save time:
do ix=0,nxm1
  do iy=0,ny
    ula(iy,ix)=pt*zzpre(iy,ix)+ptc*zz(iy,ix)
    vla(iy,ix)=pt*bbpre(iy,ix)+ptc*bb(iy,ix)
  enddo
enddo

!Write gridded fields to file:                                                                    
do ix=0,nxm1
  do iy=0,ny
    zzr4(iy,ix)=real(ula(iy,ix))
    bbr4(iy,ix)=real(vla(iy,ix))
  enddo
enddo
tr4=real(tgrid)

write(31,rec=igrids) tr4,zzr4
write(32,rec=igrids) tr4,bbr4

 !Compute domain integrals of zeta, zeta^2, b and b^2:
call average(ula,zzl1)
zzl1=zzl1*domarea
call l2norm(ula,zzl2)
call average(vla,bbl1)
bbl1=bbl1*domarea
call l2norm(vla,bbl2)

 !Write diagnostics to the files norms.dat & ene.dat:
write(13,'(f7.2,4(1x,f14.9))') tgrid,bbl1,bbl2,zzl1,zzl2
write(15,'(f7.2,4(1x,f14.9))') tgrid,ekin+apot,ekin,apot,dkdt

write(*,'(a,f7.2,4(a,f10.7))') &
    & ' t = ',tgrid,'  b_1 = ',bbl1,'  z_1 = ',zzl1,'  K = ',ekin,'  P = ',apot

 !Compute 1d vorticity spectrum:
call spec1d_fc(ula,zspec)
sumzspec=zero
do k=0,kmax
  sumzspec=sumzspec+zspec(k)
   !Normalise to take into account uneven sampling of wavenumbers 
   !in each shell [k-1/2,k+1/2]:
  zspec(k)=spmf(k)*zspec(k)
enddo
sumzspec=8.d0*sumzspec*dsumi
 !Write out spectrum to file:
write(50,'(f7.2,2(1x,f14.9),1x,i5)') tgrid,sumzspec,zzl2,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(50,'(2(1x,f12.8))') alk(k),log10(zspec(k))
enddo
 !Note: alk(k) = log_10(k)

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves zz & bb contours and residual zz for post-processing via
! congen.f90

implicit double precision(a-h,o-z)
implicit integer(i-n)

real:: zdr4(0:ny,0:nxm1),tr4
integer:: iop(max(nz,nb))
character(len=3):: pind

write(*,'(a,f12.5)') ' Saving contours at t = ',t
irec=nint(t/tcsave)
write(pind(1:3),'(i3.3)') irec

tr4=real(t)
 !Write contours to the cont subdirectory:                                                         
write(80,'(i8,1x,i9,1x,f12.5,2(1x,f16.12))') nz,nptz,t,zjump,zavg
write(90,'(i8,1x,i9,1x,f12.5,2(1x,f16.12))') nb,nptb,t,bjump,bavg

 !Save residual needed to build ultra-fine-grid vorticity with congen:
do ix=0,nxm1
  do iy=0,ny
    zdr4(iy,ix)=real(zz(iy,ix)-zc(iy,ix))
  enddo
enddo
write(83,rec=irec) tr4,zdr4

 !Save vorticity contours if any exist:                                                            
if (nz .gt. 0) then
   !First form iop; open/closed indicator:                                                         
  do j=1,nz
    iop(j)=nextz(i2z(j))/i1z(j)
     !iop = 0 for an open contour, and 1 for a closed one                                          
  enddo
  open(81,file='cont/zzindex'//pind,form='unformatted', &
      & access='direct',status='replace',recl=16*nz)
  write(81,rec=1) npz(1:nz),i1z(1:nz),indz(1:nz),iop(1:nz)
  close(81)

  open(82,file='cont/zznodes'//pind,form='unformatted', &
      & access='direct',status='replace',recl=16*nptz)
  write(82,rec=1) xz(1:nptz),yz(1:nptz)
  close(82)
endif

 !Save buoyancy contours:                                                                          
 !First form iop; open/closed indicator:                                                           
do j=1,nb
  iop(j)=nextb(i2b(j))/i1b(j)
   !iop = 0 for an open contour, and 1 for a closed one                                            
enddo
open(91,file='cont/bbindex'//pind,form='unformatted', &
    & access='direct',status='replace',recl=16*nb)
write(91,rec=1) npb(1:nb),i1b(1:nb),indb(1:nb),iop(1:nb)
close(91)

open(92,file='cont/bbnodes'//pind,form='unformatted', &
    & access='direct',status='replace',recl=16*nptb)
write(92,rec=1) xb(1:nptb),yb(1:nptb)
close(92)

return
end subroutine

!=======================================================================

 !Main end module
end module
