module evolution

! Module contains subroutines to evolve bb and zz fields according to the 
! algorithm detailed in caps.f90.

use common

implicit none

 !Vorticity field due to vorticity contours (physical)
double precision:: za(0:ny,0:nxm1)

 !Vorticity field due to vorticity contours (spectral)
double precision:: zc(0:nxm1,0:ny)

 !Dissipation operator (spectral):
double precision:: diss(0:nxm1,0:ny)

!Internal subroutine definitions (inherit global variables):

contains 

!=============================================================

subroutine advect
! Main subroutine for advecting fields and contours

implicit none

 !Local parameters:
double precision,parameter:: twistmax=2.5d0
!      twistmax: the maximum value of the time integral of |zeta|_max
!                between regularisations of the contours.
double precision,parameter:: zratmax=0.2d0
!      zratmax:  the maximum ratio r = <zd^2>/<zc^2> permitted.
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, or when r > zratmax, 
!      the code rebuilds all contours in a separate memory space.

 !Local variables:
double precision:: zdl2,zsl2,zrat
integer:: ireg

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Used for regularising contours:
twist=zero

 !Counter used for counting number of contour regularisations done:
if (nz .eq. 0) then
   !Ensure an early build of vorticity contours when vorticity is
   !initially zero:
  ireg=nregmax-1
else
  ireg=0
endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

    if (nz .gt. 0) then
       !Compute ratio of mean-square residual and contour vorticity:
      zc=zd
      call spctop_fc(nx,ny,zc,za,xfactors,yfactors,xtrig,ytrig)
      call l2norm(za,zdl2)
      call con2grid(za,xz,yz,zjump,zavg,nextz,nptz,-1)
      call l2norm(za,zsl2)
      zrat=zdl2/zsl2
    else
      zrat=zero
    endif

     !Don't continue if maximum number of regularisations reached
     !or if <zd^2>/<zc^2> > zratmax:
    if (ireg .eq. nregmax .or. zrat .gt. zratmax) then
       !Prepare residual zd for recontouring (and preserve zs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery(xb,yb,nextb,indb,npb,i1b,i2b,nb,nptb)
    call surgery(xz,yz,nextz,indz,npz,i1z,i2z,nz,nptz)
     !Record contour complexity to complexity.dat:
    write(14,'(1x,f12.5,2(1x,i8,1x,i9))') t,nb,nptb,nz,nptz

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance
  
enddo
 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Save final data if not done already:
if (int(t/tgsave) .eq. igrids) then
  call inversion
  call adapt
  call savegrid
endif
if (int(t/tcsave) .eq. iconts) call savecont

return
end subroutine advect

!=======================================================================

subroutine init
! Initialises quantities needed for normal time integration following 
! contour regeneration

implicit none

!--------------------------------------------------------------
 !Record contour complexity to complexity.dat:
write(14,'(1x,f12.5,2(1x,i8,1x,i9))') t,nb,nptb,nz,nptz

!-------------------------------------------------------------
 !Logicals used for saving gridded fields and contours:
gsave=.false.
csave=.false.

!--------------------------------------------------------------
 !Convert vorticity contours to gridded values (zc):
call con2grid(za,xz,yz,zjump,zavg,nextz,nptz,-1)
call ptospc_fc(nx,ny,za,zc,xfactors,yfactors,xtrig,ytrig)
 !zc must be in spectral space for use below; za is overwritten

 !Define (spectral) residual vorticity, zd = (1-F)[zs-zc]:
zd=bfhi*(zs-zc)
 !Here bfhi = 1-F is a high-pass spectral filter (see spectral.f90)

!--------------------------------------------------------------
 !Get the initial velocity field (uu,vv):
if (zjump .gt. zero) then 
  call inversion
else
   !Here vorticity is zero identically - no need to do inversion:
  uu=zero
  vv=zero 
  zz=zero 
endif

return
end subroutine init

!=======================================================================

subroutine prepare
! Prepares for re-contouring.  We leave this module after this routine.

implicit none

 !Local variables:
double precision:: zzmin,zzl1,zzl2
integer:: ix,iy

!--------------------------------------------------------------
 !Compute contour interval for vorticity:
zzmin=sqrt(dsumi*(f12*sum(zz(0,:)**2+zz(ny,:)**2)+sum(zz(1:nym1,:)**2)))
 !Note: zz is available from a previous call to inversion
zzl1=zero
zzl2=zero
do ix=0,nxm1
  if (abs(zz(0,ix)) .gt. zzmin) then
    zzl1=zzl1+abs(zz(0,ix))
    zzl2=zzl2+zz(0,ix)**2
  endif
  if (abs(zz(ny,ix)) .gt. zzmin) then
    zzl1=zzl1+abs(zz(ny,ix))
    zzl2=zzl2+zz(ny,ix)**2
  endif
enddo
zzl1=f12*zzl1
zzl2=f12*zzl2
do ix=0,nxm1
  do iy=1,nym1
    if (abs(zz(iy,ix)) .gt. zzmin) then
      zzl1=zzl1+abs(zz(iy,ix))
      zzl2=zzl2+zz(iy,ix)**2
    endif
  enddo
enddo
zjump=zzl2/(zzl1*dble(ncontz))
 !ncontz is set in parameters.f90. 

!--------------------------------------------------------------
 !Convert vorticity contours to gridded values (zc):
call con2grid(za,xz,yz,zjump,zavg,nextz,nptz,-1)
call ptospc_fc(nx,ny,za,zc,xfactors,yfactors,xtrig,ytrig)
 !zc must be in spectral space for use below; za is overwritten

 !Combine spectral vorticity anomaly (zs) and residual (zd):
zs=bflo*zs+bfhi*zc+zd
zd=zs-zc
 !Average values are not used here:
zs(0,0)=zero
zd(0,0)=zero

 !Convert zd to spectral space as zz for recontouring:
call spctop_fc(nx,ny,zd,zz,xfactors,yfactors,xtrig,ytrig)
 !zd is overwritten, but it is redefined in init upon return to this module
 !zz is also re-computed in init upon return.

return
end subroutine prepare

!=======================================================================

subroutine inversion
! Inverts Laplace's operator on vorticity (zz) to obtain the 
! streamfunction (pp) and the velocity (uu,vv) = (-dpp/dy,dpp/dx).

implicit none

 !Local variable (physical):
double precision:: za(0:ny,0:nxm1)

!------------------------------------------------------------
 !Call con2grid to get updated contour vorticity (zc):
call con2grid(za,xz,yz,zjump,zavg,nextz,nptz,-1)
call ptospc_fc(nx,ny,za,zc,xfactors,yfactors,xtrig,ytrig)
 !zc must be in spectral space for use below; za is overwritten

 !Combine vorticity fields to update zc with full (spectral) field,
 !F[zs-zc]+zc+zd -> zc, where F = bflo is a low pass filter:
zc=bflo*zs+bfhi*zc+zd !bfhi = 1 - F here (see spectral.f90)
zc(0,0)=zero !Average value is set in main_invert below

 !Invert Laplace operator on zc spectrally and obtain velocity field: 
call main_invert(zc,uu,vv,zz,zavg)
 !Main invert returns also the gridded vorticity zz, and corrects its
 !average to be zavg.

return
end subroutine inversion

!=======================================================================

subroutine advance

! Advances fields from time t to t+dt using an iterative implicit 
! trapezoidal method of the form
!
!     (F^{n+1}-F^n)/dt = (L^{n+1}+L^n)/2 + (S^{n+1}+S^n)/2
!
! for a field F, where n refers to the time level, L[F] refers to
! the linear dissipation terms (hyperdiffusion), and S[F] refers to
! the remaining source terms.

! We start with the guess S^{n+1} = S^n and iterate  niter  times
! (see parameter statement below).

! Contours are treated similarly except L = 0 and S = (u,v).

implicit none

 !Number of iterations of above scheme:
integer,parameter:: niter=2

 !Spectral fields:
double precision:: zsi(0:nxm1,0:ny),szs(0:nxm1,0:ny)
double precision:: zdi(0:nxm1,0:ny),zdm(0:nxm1,0:ny),szd(0:nxm1,0:ny)

 !Contour arrays:
double precision:: ub(nptb),vb(nptb),xbm(nptb),ybm(nptb)
double precision:: uz(nptz),vz(nptz),xzm(nptz),yzm(nptz)

 !Other local quantities:
double precision:: xx
integer:: i,iter

!-------------------------------------------------------------------
 !Invert vorticity for velocity at current time level, say t=t^n:
call inversion

 !Adapt the time step and save various diagnostics each time step:
call adapt

 !Possibly save data (gsave & csave set by adapt):
if (gsave) call savegrid
if (csave) call savecont

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all fields:

 !Calculate the vorticity source terms:
call zztend(szd,szs)

 !Initialise iteration (dt2 = dt/2 and dt4 = dt/4 below):
zsi=zs+dt2*szs
zs=zs+dt*szs
zdi=zd
zdm=zd+dt4*szd
zd=diss*(zdm+dt4*szd)-zdi
 !diss is related to the hyperdiffusive operator (see end of adapt)

 !Contours:
call velint(uu,vv,xb,yb,ub,vb,nptb)
do i=1,nptb
  xx=xb(i)+dt2*ub(i)
  xbm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  ybm(i)=min(ymax,max(ymin,yb(i)+dt2*vb(i)))
  xx=xb(i)+dt*ub(i)
  xb(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yb(i)=min(ymax,max(ymin,yb(i)+dt*vb(i)))
enddo

if (nptz .gt. 0) then
  call velint(uu,vv,xz,yz,uz,vz,nptz)
  do i=1,nptz
    xx=xz(i)+dt2*uz(i)
    xzm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yzm(i)=min(ymax,max(ymin,yz(i)+dt2*vz(i)))
    xx=xz(i)+dt*uz(i)
    xz(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yz(i)=min(ymax,max(ymin,yz(i)+dt*vz(i)))
  enddo
endif

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion

   !Calculate the vorticity source terms:
  call zztend(szd,szs)

   !Update fields:
  zs=zsi+dt2*szs
  zd=diss*(zdm+dt4*szd)-zdi

   !Update contours:
  call velint(uu,vv,xb,yb,ub,vb,nptb)
  do i=1,nptb
    xx=xbm(i)+dt2*ub(i)
    xb(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yb(i)=min(ymax,max(ymin,ybm(i)+dt2*vb(i)))
  enddo

  if (nptz .gt. 0) then
    call velint(uu,vv,xz,yz,uz,vz,nptz)
    do i=1,nptz
      xx=xzm(i)+dt2*uz(i)
      xz(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
      yz(i)=min(ymax,max(ymin,yzm(i)+dt2*vz(i)))
    enddo
  endif
enddo

 !Update gridded buoyancy field (bb) from contours (xb,yb):
call con2grid(bb,xb,yb,bjump,bavg,nextb,nptb,-1)

 !Advance time:
t=t+dt

return
end subroutine advance

!=======================================================================

subroutine zztend(szd,szs)
! Computes the residual and full vorticity local rates-of-change
! in spectral space

implicit none

 !Passed variables:
double precision:: szd(0:nxm1,0:ny),szs(0:nxm1,0:ny)

 !Local variables (physical):
double precision:: px(0:ny,0:nxm1),py(ny,0:nxm1)

 !Local variables (spectral):
double precision:: sx(0:nxm1,0:ny),sy(0:nxm1,ny)

 !Others:
double precision:: dzdt(0:ny,0:nxm1)
double precision:: wks(0:nxm1,0:ny)
double precision:: dzdxs(0:nxm1,0:ny),dzdxp(0:ny,0:nxm1)
double precision:: dzdys(0:nxm1,ny),dzdyp(ny,0:nxm1)

!--------------------------------------------------------------
 !Residual vorticity source zd_t = bb_x - (u,v)*grad(zd):

 !Get bb_x directly from buoyancy contours (see contours.f90):
call getzzsrc(dzdt) !dzdt is in physical space

 !Obtain x & y derivatives of zd -> (px,py) in physical space:
call xderiv_fc(nx,ny,hrkx,zd,sx)
call spctop_fc(nx,ny,sx,px,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,zd,sy)
call spctop_fs(nx,ny,sy,py,xfactors,yfactors,xtrig,ytrig)

 !Compute (u,v)*grad(zd) -> px in physical space:
px(0,:)=uu(0,:)*px(0,:)
px(1:nym1,:)=uu(1:nym1,:)*px(1:nym1,:)+vv(1:nym1,:)*py(1:nym1,:)
px(ny,:)=uu(ny,:)*px(ny,:)

 !Combine bb_x and (u,v)*grad(zd) to form zd_t:
px=dzdt-px

 !Convert to spectral space as szd, then apply de-aliasing filter:
call ptospc_fc(nx,ny,px,szd,xfactors,yfactors,xtrig,ytrig)
szd=filt*szd

!--------------------------------------------------------------
 !Conserved vorticity source zs_t = - (u,v)*grad(zs):

 !Obtain x & y derivatives of zs -> (px,py) in physical space:
call xderiv_fc(nx,ny,hrkx,zs,sx)
call spctop_fc(nx,ny,sx,px,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,zs,sy)
call spctop_fs(nx,ny,sy,py,xfactors,yfactors,xtrig,ytrig)

 !Compute (u,v)*grad(zs) -> px in physical space:
px(0,:)=uu(0,:)*px(0,:)
px(1:nym1,:)=uu(1:nym1,:)*px(1:nym1,:)+vv(1:nym1,:)*py(1:nym1,:)
px(ny,:)=uu(ny,:)*px(ny,:)

 !Convert to spectral space as szs, then apply de-aliasing filter:
call ptospc_fc(nx,ny,px,szs,xfactors,yfactors,xtrig,ytrig)
szs=-filt*szs

return
end subroutine zztend

!=======================================================================

subroutine adapt
! Adapts the time step to ensure dt < dtfac/max(|zeta|_max)

implicit none

 !For defining the max strain & buoyancy frequency based time step:
double precision,parameter:: alpha=0.1d0
 !Note: EPIC-2D paper recommends alpha = 0.2 for ls-rk4 method

 !For controlling numerical stability (CFL_max <= 0.8 recommended):
double precision,parameter:: cflmax=0.8d0
double precision,parameter:: cflpf=cflmax*glmin

 !Local variables (physical):
double precision:: px(0:ny,0:nxm1),py(ny,0:nxm1)

 !Local variables (spectral):
double precision:: sx(0:nxm1,0:ny),sy(0:nxm1,ny)

 !Other variables:
double precision:: bfmax,zzmax,zzrms,zztmp,zzl1,zzl2,zzch
double precision:: ggmax,uumax,dfac
integer:: ix,iy

!----------------------------------------------------------
 !Obtain x & y derivatives of buoyancy bb -> px, py (physical):
za=bb !Copy bb to avoid overwriting it in FFT
call ptospc_fc(nx,ny,za,zc,xfactors,yfactors,xtrig,ytrig)
call xderiv_fc(nx,ny,hrkx,zc,sx) !zc = spectral bb here
call spctop_fc(nx,ny,sx,px,xfactors,yfactors,xtrig,ytrig) !px = bb_x
call yderiv_fc(nx,ny,rky,zc,sy)
call spctop_fs(nx,ny,sy,py,xfactors,yfactors,xtrig,ytrig) !py = bb_y

 !Compute |grad{bb}|^2 = px^2 + py^2 -> px in physical space:
px(0,:)=px(0,:)**2
px(1:nym1,:)=px(1:nym1,:)**2+py(1:nym1,:)**2
px(ny,:)=px(ny,:)**2

 !Maximum buoyancy frequency:
bfmax=sqrt(sqrt(maxval(px)))

 !Maximum vorticity:
px=zz**2
zzmax=sqrt(maxval(px))

 !R.m.s. vorticity:
zzrms=sqrt(dsumi*(f12*sum(px(0,:)+px(ny,:))+sum(px(1:nym1,:))))

 !Characteristic vorticity, <zz^2>/<|zz|> for |zz| > zz_rms:
zzl1=small
zzl2=zero
do ix=0,nxm1
  do iy=1,ny
    zztmp=f12*(zz(iy-1,ix)+zz(iy,ix))
    if (abs(zztmp) .gt. zzrms) then
      zzl1=zzl1+abs(zztmp)
      zzl2=zzl2+zztmp**2
    endif
  enddo
enddo
zzch=zzl2/zzl1

 !Compute x derivative of velocity components:
za=uu
call ptospc_fc(nx,ny,za,zc,xfactors,yfactors,xtrig,ytrig)
call xderiv_fc(nx,ny,hrkx,zc,sx)
call spctop_fc(nx,ny,sx,za,xfactors,yfactors,xtrig,ytrig)
 !za = u_x
px=vv
call ptospc_fc(nx,ny,px,zc,xfactors,yfactors,xtrig,ytrig)
call xderiv_fc(nx,ny,hrkx,zc,sx)
call spctop_fc(nx,ny,sx,px,xfactors,yfactors,xtrig,ytrig)
 !px = v_x

 !Strain rate squared, u_x^2 + (v_x - zz/2)^2:
za=za**2+(px-f12*zz)**2

 !Maximum strain rate:
ggmax=sqrt(maxval(za))

 !Maximum speed:
uumax=sqrt(maxval(uu**2+vv**2))

 !Choose new time step:
dt=min(alpha/(ggmax+small),alpha/(bfmax+small),cflpf/(uumax+small),tgsave/four)

 !Update fractional time steps:
dt2=dt/two
dt4=dt/four

!---------------------------------------------------------------------
 !Update hyperdiffusion operator used in time stepping:
dfac=zzch*dt2
 !zzch is the characteristic vorticity defined above.
diss=two/(one+dfac*hdis)
 !hdis = C*(K/K_max)^{2p} where K^2 = k_x^2+k_y^2, p is the order,
 !K_max is the maximum x or y wavenumber and C is a dimensionless
 !prefactor (see spectral.f90 and parameters.f90 where C = prediss).

 !Increment the integral of |zeta|_max:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
 !Save |u|_max, N_max and gamma_max to monitor.asc:
write(22,'(1x,f13.6,3(1x,1p,e14.7))') t,uumax,bfmax,ggmax

 !Save vorticity diagnostics to vorticity.asc:
write(23,'(1x,f13.6,3(1x,1p,e14.7))') t,zzmax,zzrms,zzch

!---------------------------------------------------------------------
 !Set flags to save data:
gsave=(int(t/tgsave) .eq. igrids)
 !Gridded data will be saved at time t if gsave is true.
csave=(int(t/tcsave) .eq. iconts)
 !Contour data will be saved at time t if csave is true.

return
end subroutine adapt

!=======================================================================
      
subroutine savegrid
! Saves vorticity and buoyancy fields, energy components, as well as
! 1d spectra.

implicit none

 !Local variables:
double precision:: spec(0:max(nx,ny))
double precision:: eke
integer:: k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Compute available potential energy (ape) by calling con2grid:
call con2grid(bb,xb,yb,bjump,bavg,nextb,nptb,iene)
 !If iene = 0, this also computes the reference PE (usually at t = 0).
iene=1

 !Compute kinetic energy
za=uu**2+vv**2
eke=f12*garea*(f12*sum(za(0,:)+za(ny,:))+sum(za(1:nym1,:)))

 !Write energy:
write(21,'(1x,f13.6,3(1x,1p,e14.7))') t,eke,ape,eke+ape

write(*,'(a,f7.2,3(a,1p,e14.7))') &
    & ' t = ',t,'  K = ',eke,'  P = ',ape,'  E = ',eke+ape

 !Write gridded vorticity and buoyancy:
write(31,rec=igrids) real(t),real(zz)
write(32,rec=igrids) real(t),real(bb)

!---------------------------------------------------------------
 !Compute 1d vorticity & buoyancy spectra:

za=zz !Copy zz to avoid overwriting it in FFT
call ptospc_fc(nx,ny,za,zc,xfactors,yfactors,xtrig,ytrig)
call spec1d_fc(zc,spec) !zc is used temporarily here
spec(0:kmax)=log10(spmf(0:kmax)*spec(0:kmax)+1.d-32)
write(51,'(f13.6,1x,i5)') t,kmax
do k=1,kmax
  write(51,'(3(1x,f12.8))') alk(k),spec(k)
enddo

 !Compute 1d buoyancy spectrum and write:
za=bb !Copy bb to avoid overwriting it in FFT
call ptospc_fc(nx,ny,za,zc,xfactors,yfactors,xtrig,ytrig)
call spec1d_fc(zc,spec) !zc is used temporarily here
spec(0:kmax)=log10(spmf(0:kmax)*spec(0:kmax)+1.d-32)
write(52,'(f13.6,1x,i5)') t,kmax
do k=1,kmax
  write(52,'(3(1x,f12.8))') alk(k),spec(k)
enddo

 !spmf takes into account uneven sampling of wavenumbers in each
 !shell [K-1/2,K+1/2], where K is proportional to k (see spectral.f90)

 !alk(k) = log_10(K)

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont
! Saves zz & bb contours and residual zz for post-processing

implicit none

 !Local variables:
integer:: iop(max(nz,nb)),j
character(len=3):: pind

!---------------------------------------------------------------
 !Increment counter for direct file access:
iconts=iconts+1

write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') iconts-1

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,2(1x,f16.12))') nz,nptz,t,zjump,zavg
write(90,'(i8,1x,i9,1x,f12.5,2(1x,f16.12))') nb,nptb,t,bjump,bavg

 !Save residual needed to build ultra-fine-grid vorticity with congen:
zc=zd !Copy zd to avoid overwriting it
call spctop_fc(nx,ny,zc,za,xfactors,yfactors,xtrig,ytrig)
write(83,rec=iconts) real(t),real(za)

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

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
