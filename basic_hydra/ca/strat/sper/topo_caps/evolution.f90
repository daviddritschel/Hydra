module evolution

! Module contains subroutines to evolve bb and zz fields according to the 
! algorithm detailed in casl.f90.

use common

implicit none

 !Energies:
double precision:: ekepre,ekepost,epepre,epepost

 !Physical fields:
double precision:: zz(0:ny,0:nxm1),zc(0:ny,0:nxm1),szs(0:ny,0:nxm1),szd(0:ny,0:nxm1)
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: zzpre(0:ny,0:nxm1),zspre(0:ny,0:nxm1),zdpre(0:ny,0:nxm1)
double precision:: bbpre(0:ny,0:nxm1),szspre(0:ny,0:nxm1),szdpre(0:ny,0:nxm1)

 !Spectral arrays:
double precision:: zopi(0:nxm1,0:ny)

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
    call con2grid(zc,xz,yz,zjump,nextz,nptz)
     !Reset zs and zd:
    call reset(zc,zs,zd,zavg)

     !Copy gridded vorticity fields to old time level:
    do ix=0,nxm1
      do iy=0,ny
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
      szspre(iy,ix)=szs(iy,ix)
      szdpre(iy,ix)=szd(iy,ix)
      zdpre(iy,ix)=zd(iy,ix)
      zspre(iy,ix)=zs(iy,ix)
      zzpre(iy,ix)=zz(iy,ix) 
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

 !Work arrays:
double precision:: wka(0:ny,0:nxm1),wkb(0:ny,0:nxm1)

!---------------------------------------------
 !Record contour complexity to complexity.dat:
write(14,'(1x,f12.5,2(1x,i8,1x,i9))') t,nb,nptb,nz,nptz

!---------------------------------------------------
 !Convert vorticity contours to gridded values (zc):
call con2grid(zc,xz,yz,zjump,nextz,nptz)

 !Define residual vorticity zd = zs-zc-F[zs-zc]
do ix=0,nxm1
  do iy=0,ny
    wka(iy,ix)=zs(iy,ix)-zc(iy,ix)
    wkb(iy,ix)=wka(iy,ix)
  enddo
enddo

call filter(wka,0,2)

do ix=0,nxm1
  do iy=0,ny
    zd(iy,ix)=wkb(iy,ix)-wka(iy,ix)
  enddo
enddo
 !Ensure zero domain average zd at beginning of next time step:
call restore(zd,zero)

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

 !Calculate the source terms for zd and zs at time t:
call zztend(szd,szs,t)

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
call con2grid(zc,xz,yz,zjump,nextz,nptz)

 !Combine fields to update zz with full field:
call combine(zz,zc,zs,zd,zavg)

 !Invert Laplace operator on zz spectrally and obtain velocity field: 
call main_invert(zz,uu,vv,t)      

return
end subroutine

!=======================================================================
      
subroutine advance

! Computes bb(t+dt) by contour advection and zz(t+dt) by a combination of 
! contour advection and a standard pseudo-spectral scheme (PS)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
double precision:: ub(nptb),vb(nptb),xbm(nptb),ybm(nptb)
double precision:: uz(nptz),vz(nptz),xzm(nptz),yzm(nptz)

!------------------------------------------------------------------------
 !Calculate the source term for zeta at time t:
call zztend(szd,szs,t)

 !Copy gridded vorticity & form vorticity tendency term at old time level (time t):
do ix=0,nxm1
  do iy=0,ny
    zdpre(iy,ix)=zd(iy,ix) 
    zspre(iy,ix)=zs(iy,ix) 
    szdpre(iy,ix)=zd(iy,ix)+qudt*szd(iy,ix)
    szspre(iy,ix)=zs(iy,ix)+qudt*szs(iy,ix)
  enddo
enddo
 !At the end of the subroutine zz will contain vorticity at t + dt

 !Prepare contour evolution; get velocity on buoyancy contour nodes:
call velint(uu,vv,xb,yb,ub,vb,nptb)
do i=1,nptb
  xx=xb(i)+hfdt*ub(i)
  xbm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  ybm(i)=min(ymax,max(ymin,yb(i)+hfdt*vb(i)))
  xx=xb(i)+dt*ub(i)
  xb(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yb(i) =min(ymax,max(ymin,yb(i)+dt*vb(i)))
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
  call ps_step

   !Obtain zz & hence uu & vv at time t + dt:
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

   !Update vorticity tendency:
  call zztend(szd,szs,t+dt)
enddo

 !Obtain final corrected zs & zd at time t + dt:
call ps_step

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
call con2grid(zc,xz,yz,zjump,nextz,nptz)

 !Combine fields to update zz with full field:
call combine(zz,zc,zs,zd,zavg)

return
end subroutine

!=======================================================================
      
subroutine ps_step

! Evolve zd & zs from t to t + dt using pseudo-spectral method.

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: zdbs(0:nxm1,0:ny)

do ix=0,nxm1
  do iy=0,ny
    zd(iy,ix)=szdpre(iy,ix)+qudt*szd(iy,ix)
  enddo
enddo

call ptospc_fc(nx,ny,zd,zdbs,xfactors,yfactors,xtrig,ytrig)

do ky=0,ny
  do kx=0,nxm1
    zdbs(kx,ky)=zopi(kx,ky)*zdbs(kx,ky)
  enddo
enddo

call spctop_fc(nx,ny,zdbs,zd,xfactors,yfactors,xtrig,ytrig)

do ix=0,nxm1
  do iy=0,ny
    zd(iy,ix)=two*zd(iy,ix)-zdpre(iy,ix)
    zs(iy,ix)=two*(szspre(iy,ix)+qudt*szs(iy,ix))-zspre(iy,ix)
  enddo
enddo

return
end subroutine

!=======================================================================
      
subroutine zztend(vard,vars,time)

! Get vorticity local rate-of-change \pa{zs}/\pa{t} & \pa{zd}/\pa{t}:

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: vard(0:ny,0:nxm1),vars(0:ny,0:nxm1),ztmp(0:ny,0:nxm1)
double precision:: wks(0:nxm1,0:ny)
double precision:: dzdxs(0:nxm1,0:ny),dzdxp(0:ny,0:nxm1)
double precision:: dzdys(0:nxm1,ny),dzdyp(ny,0:nxm1)

!--------------------------------------------------------------
 !Get linear source term (Dzz/Dt):
call getzzsrc(vard,time,dyoridx,dyoridy)

 !copy zd to avoid overwriting:
do ix=0,nxm1
  do iy=0,ny
    ztmp(iy,ix)=zd(iy,ix)
  enddo
enddo

 !Obtain x & y derivatives of vorticity:
call ptospc_fc(nx,ny,ztmp,wks,xfactors,yfactors,xtrig,ytrig)

 !Apply de-aliasing filter:
do ky=0,ny
  do kx=0,nxm1
    wks(kx,ky)=wks(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo

call xderiv_fc(nx,ny,hrkx,wks,dzdxs)
call spctop_fc(nx,ny,dzdxs,dzdxp,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,wks,dzdys)
call spctop_fs(nx,ny,dzdys,dzdyp,xfactors,yfactors,xtrig,ytrig)

 !Subtract u.grad(zd) from the linear source term to get \pa{zd}/\pa{t}:
do ix=0,nxm1
  do iy=1,nym1
    vard(iy,ix)=vard(iy,ix)-uu(iy,ix)*dzdxp(iy,ix)-vv(iy,ix)*dzdyp(iy,ix)
  enddo
enddo
 !Deal with edge values separately:
do ix=0,nxm1
  vard(0, ix)=vard(0 ,ix)-uu(0 ,ix)*dzdxp(0 ,ix)
  vard(ny,ix)=vard(ny,ix)-uu(ny,ix)*dzdxp(ny,ix)
enddo

!------------------------------------------------------
 !copy zs to avoid overwriting:
do ix=0,nxm1
  do iy=0,ny
    ztmp(iy,ix)=zs(iy,ix)
  enddo
enddo

 !Obtain x & y derivatives of vorticity:
call ptospc_fc(nx,ny,ztmp,wks,xfactors,yfactors,xtrig,ytrig)

 !Apply de-aliasing filter:
do ky=0,ny
  do kx=0,nxm1
    wks(kx,ky)=wks(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo

call xderiv_fc(nx,ny,hrkx,wks,dzdxs)
call spctop_fc(nx,ny,dzdxs,dzdxp,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,wks,dzdys)
call spctop_fs(nx,ny,dzdys,dzdyp,xfactors,yfactors,xtrig,ytrig)

 !Subtract u.grad(zs) from the linear source term to get \pa{zs}/\pa{t}:
do ix=0,nxm1
  do iy=1,nym1
    vars(iy,ix)=-uu(iy,ix)*dzdxp(iy,ix)-vv(iy,ix)*dzdyp(iy,ix)
  enddo
enddo
 !Deal with edge values separately:
do ix=0,nxm1
  vars(0, ix)=-uu(0 ,ix)*dzdxp(0 ,ix)
  vars(ny,ix)=-uu(ny,ix)*dzdxp(ny,ix)
enddo

return
end subroutine

!=======================================================================

subroutine adapt(igsave,icsave)

! Adapts the time step to ensure dt < dtfac/max(|zeta|_max)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Work arrays:
double precision:: wka(0:ny,0:nx),wkb(0:ny,0:nx)

 !Maximum CFL number:
double precision,parameter:: cflmax=0.7d0,glmin=min(glx,gly)

 !Compute accurate advection time step (dtacc):
zzmax=small
uumax=small
do ix=0,nxm1
  do iy=0,ny
    zzmax=max(zzmax,abs(zz(iy,ix)))
    uumax=max(uumax,uu(iy,ix)**2+vv(iy,ix)**2)
  enddo
enddo
dtacc=dtfac/zzmax
uumax=sqrt(uumax)
dtcfl=cflmax*glmin/uumax

!---------------------------------------------------------------------
 !Choose a new time step: 
if (dt .gt. zero) then 
  dt=min(dtacc,dtcfl,dtmax)
  if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
else
   !Calculate a maximum timestep (dtmax) based on grad(bb):
  call bxderiv(bb,wka)
  call byderiv(bb,wkb)
  dbmax=small
  do ix=0,nxm1
    do iy=0,ny
      dbmax=max(dbmax,confaci(iy,ix)*(wka(iy,ix)**2+wkb(iy,ix)**2))
    enddo
  enddo 
   !Limit max timestep to a data save time
  dtmax=min(tgsave,tcsave,dtfac/sqrt(sqrt(dbmax)))
  dt=min(dtacc,dtcfl,dtmax)
endif 
hfdt=dt/two
qudt=dt/four

!---------------------------------------------------------------------
 !Define spectral evolution operator (ignoring conformal factor):
if (nnu .eq. 1) then   
   !Laplacian dissipation
  do ky=0,ny
    do kx=0,nxm1
      zopi(kx,ky)=one/(one+hfdt*diss(kx,ky))
    enddo
  enddo
else
   !Hyper-viscous dissipation (proportional to |zz|_max)
  dissfac=hfdt*zzmax   
  do ky=0,ny
    do kx=0,nxm1
      zopi(kx,ky)=one/(one+dissfac*diss(kx,ky))
    enddo
  enddo
endif

!---------------------------------------------------------------------
 !Compute CFL number:
cfl=uumax*dt/glmin

 !Record cfl, max |z|, and max and l1,l2 norms of b to monitor.dat:
write(12,'(1x,f12.5,1x,f6.4,1x,1p,e14.7)') t,cfl,zzmax

!---------------------------------------------------------------------

 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
tgrid=tgsave*dble(itime)+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  igsave=1
   !Compute current gridded buoyancy and store in bbpre:
  call con2grid(bbpre,xb,yb,bjump,nextb,nptb)
   !Restore correct average:
  call restore(bbpre,bavg)
   !Compute kinetic energy:
  call kinetic(uu,vv,ekepre)
   !Compute potential energy:
  call potential(bbpre,epepre)
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

real:: zzr4(0:ny,0:nxm1),bbr4(0:ny,0:nxm1),tr4

 !Work arrays:
double precision:: wka(0:ny,0:nxm1),wkb(0:ny,0:nxm1)

 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Update gridded buoyancy field:
call con2grid(bb,xb,yb,bjump,nextb,nptb)

 !Restore correct average:
call restore(bb,bavg)

 !Compute kinetic energy
call kinetic(uu,vv,ekepost)

 !Compute potential energy
call potential(bb,epepost)

 !Compute time interpolated energies:
ekin=pt*ekepre+ptc*ekepost
epot=pt*epepre+ptc*epepost

 !Store vorticity and buoyancy at save time:
do ix=0,nxm1
  do iy=0,ny
    wka(iy,ix)=pt*zzpre(iy,ix)+ptc*zz(iy,ix)
    wkb(iy,ix)=pt*bbpre(iy,ix)+ptc*bb(iy,ix)
  enddo
enddo

!Write gridded fields to file:                                                                    
do ix=0,nxm1
  do iy=0,ny
    zzr4(iy,ix)=real(wka(iy,ix))
    bbr4(iy,ix)=real(wkb(iy,ix))
  enddo
enddo
tr4=real(tgrid)

write(31,rec=igrids) tr4,zzr4
write(32,rec=igrids) tr4,bbr4

 !Compute domain integrals of zeta, zeta^2, b and b^2:
call average(wka,zzl1)
zzl1=zzl1*domarea
call l2norm(wka,zzl2)
call average(wkb,bbl1)
bbl1=bbl1*domarea
call l2norm(wkb,bbl2)

 !Write diagnostics to the files norms.dat & ene.dat:
write(13,'(f7.2,4(1x,f14.9))') tgrid,bbl1,bbl2,zzl1,zzl2
write(15,'(f7.2,3(1x,f14.9))') tgrid,ekin,epot,ekin+epot

write(*,'(a,f7.2,4(a,f9.6))') ' t = ',tgrid, &
     '  b_1 = ',bbl1,'  z_1 = ',zzl1,'  K = ',ekin,'  APE = ',epot

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves zz & bb contours and residual zz for post-processing via
! congen.f90

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: wka(0:ny,0:nxm1)
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
    wka(iy,ix)=zz(iy,ix)-zc(iy,ix)
  enddo
enddo
 !Ensure zero average:
call restore(wka,zero)
 !Convert to real*4 and write:
zdr4=real(wka)
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
