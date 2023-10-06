module evolution

! Module contains subroutines to evolve PV contours and all fields 
! according to the algorithm detailed in casl.f90.

use common

implicit none

 !Energies:
double precision:: eneupre,eneupost,enebpre,enebpost

 !Physical fields:
double precision:: qq(ny,nx),qc(ny,nx),sq(ny,nx)
double precision:: qspre(ny,nx),qdpre(ny,nx),sqpre(ny,nx)
double precision:: uu(ny,nx),vv(ny,nx),pp(ny,nx)
double precision:: aapre(ny,nx),aax(ny,nx),aay(ny,nx)
 !Semi-Lagrangian advection arrays:
double precision:: x0(ny,nx),y0(ny,nx),ula(ny,nx),vla(ny,nx)
 !Allocatable source arrays:
double precision,allocatable,dimension(:,:):: sqstoch

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

integer:: ireg,igsave,icsave,ix,iy

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Used for regularising contours:
twist=zero

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

   !Perform surgery & field reset when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1
     !Don't continue if maximum number of regularisations reached;
     !it is time to recontour (pass control back to main program):
    if (ireg .eq. nregmax) then
       !Prepare PV residual for recontouring:
      do ix=1,nx
        do iy=1,ny
          qd(iy,ix)=qq(iy,ix)-qc(iy,ix)
          qs(iy,ix)=qq(iy,ix)
        enddo
      enddo
       !Re-compute PV contour interval (done only if beta = 0):
      if (.not. beffect) call contint(qq,ncontq,qjump)
       !De-allocate memory for the stochastic PV source array if used:
      if (stoch) deallocate(sqstoch)
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

     !Convert PV contours to gridded values (qc):
    call con2grid(qc)

    twist=twist-twistmax
  endif

   !Adjust timestep (dt) on maximum vorticity magnitude:
  call adapt(igsave,icsave)

   !Advect flow from time t to t + dt:
  call advance

   !Update the time:
  t=t+dt

   !Possibly save fields & energy at chosen save time (tgrid):
  if (igsave .eq. 1) call savegrid

   !Possibly save contours for post processing:
  if (icsave .eq. 1) call savecont
  
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
 !Record contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

 !Allocate memory for stochastic PV source array if used:
if (stoch) allocate(sqstoch(ny,nx))

 !Convert PV contours to gridded values (qc):
call con2grid(qc)

 !Define residual PV qd = qs-qc-F[qs-qc]:
do ix=1,nx
  do iy=1,ny
    ula(iy,ix)=qs(iy,ix)-qc(iy,ix)
    vla(iy,ix)=ula(iy,ix)
  enddo
enddo

call filter(ula,2)

do ix=1,nx
  do iy=1,ny
    qd(iy,ix)=vla(iy,ix)-ula(iy,ix)
    qq(iy,ix)=qs(iy,ix)
  enddo
enddo

 !Get the initial velocity field (uu,vv):
if (qjump .gt. zero) then 
  call main_invert(qq,uu,vv,pp,.true.)
else
   !Here PV is zero identically - no need to do inversion:
  do ix=1,nx
    do iy=1,ny
      uu(iy,ix)=zero
      vv(iy,ix)=zero 
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine advance

! Computes qq(t+dt) by contour advection by trajectory 
! integration using a standard semi-Lagrangian (SL) scheme.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
double precision:: aad(ny,nx)
double precision:: uq(nptq),vq(nptq),xqm(nptq),yqm(nptq)

!-------------------------------------------------------------------
 !Reset qs and qd at beginning of each time step (qc assumed known):
call reset(qc,qs,qd)

 !Calculate the source term for PV:
call source(sq)

 !Copy gridded semi-Lagrangian fields to old time level:
do ix=1,nx
  do iy=1,ny
    qdpre(iy,ix)=qd(iy,ix) 
    qspre(iy,ix)=qs(iy,ix) 
    aapre(iy,ix)=aa(iy,ix)
    sqpre(iy,ix)=sq(iy,ix)
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
do ix=1,nx
  do iy=1,ny
    x0(iy,ix)=mod(xigmax+xig(ix)-gcx*uu(iy,ix),xigmax)
    y0(iy,ix)=mod(yigmax+yig(iy)-gcy*vv(iy,ix),yigmax)
  enddo
enddo
 !Note, (uu,vv) is used since we have no other velocity field available

 !Compute a slightly diffused form of A used in sl_step below:
call diffuse(aa,aad,hfdt)

 !Prepare contour evolution; get velocity on PV contour nodes:
if (nptq .gt. 0) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xq(i)+hfdt*uq(i)
    yy=yq(i)+hfdt*vq(i)
    xqm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqm(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    xx=xq(i)+dt*uq(i)
    yy=yq(i)+dt*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

 !Iterate to converge on implicit trapezoidal integration:
do iter=1,niter
   !Obtain aa, qs & qd at time t + dt:
  call sl_step(aad)

   !Obtain qq & hence uu & vv at time t + dt:
  call inversion(damping)

   !Correct departure grid point (x0,y0):
  do ix=1,nx
    do iy=1,ny
       !Obtain old velocity (time t) at the departure grid point using
       !bi-linear interpolation of (ula,vla):
      ix0=1+int(x0(iy,ix))
      pxc=dble(ix0)-x0(iy,ix)
      px=one-pxc
      ix1=ixp(ix0)

      iy0=1+int(y0(iy,ix))
      pyc=dble(iy0)-y0(iy,ix)
      py=one-pyc
      iy1=iyp(iy0)

      uod=pyc*(pxc*ula(iy0,ix0)+px*ula(iy0,ix1)) &
      &   +py*(pxc*ula(iy1,ix0)+px*ula(iy1,ix1))

      vod=pyc*(pxc*vla(iy0,ix0)+px*vla(iy0,ix1)) &
      &   +py*(pxc*vla(iy1,ix0)+px*vla(iy1,ix1))

      x0(iy,ix)=mod(xigmax+xig(ix)-hgcx*(uod+uu(iy,ix)),xigmax)
      y0(iy,ix)=mod(yigmax+yig(iy)-hgcy*(vod+vv(iy,ix)),yigmax)
       !(uu,vv) is the new velocity (time t+dt) at the arrival grid point
    enddo
  enddo

   !Update the PV contour points:
  if (nptq .gt. 0) then
    call velint(uu,vv,uq,vq)
    do i=1,nptq
      xx=xqm(i)+hfdt*uq(i)
      yy=yqm(i)+hfdt*vq(i)
      xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
      yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    enddo
  endif

   !Update source terms for PV using latest fields at t + dt:
  call source(sq)

enddo

 !Obtain final corrected aa, qs & qd at time t + dt:
call sl_step(aad)

 !Obtain final corrected uu & vv at time t + dt from qq:
call inversion(.true.)

 !Update the PV contour points:
if (nptq .gt. 0) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqm(i)+hfdt*uq(i)
    yy=yqm(i)+hfdt*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

 !Call con2grid to get updated contour PV (qc):
call con2grid(qc)

 !Combine fields to update qq with full field:
call combine(qq,qc,qs,qd)

return
end subroutine

!=======================================================================

subroutine sl_step(aad)

! Carries out semi-Lagrangian integration of aa, qs & qd from t to t+dt

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed array:
double precision:: aad(ny,nx)
 !Local arrays:
double precision:: ddy(ny,nx),faa(ny,nx)

!----------------------------------------------------------------------
 !Integrate qs using bi-cubic Lagrange interpolation of qspre at x0,y0:
call interpol(qspre,qs,x0,y0)
 !This approximates Dqs/Dt = 0 when beta = 0

 !Integrate A carrying forward the diffused field A_d:
call interpol(aad,aa,x0,y0)

 !Add B_0*(y-yg) to incorporate linear background gradient in total A:
if (beffect) then
   !Add also displacement in y to qs to incorporate beta*y:
  do ix=1,nx
    do iy=1,ny
      yy=gly*(y0(iy,ix)-yig(iy))
      ddy(iy,ix)=yy-elly*dble(int(yy*hlyi))
      qs(iy,ix)=qs(iy,ix)+beta*ddy(iy,ix)
    enddo
  enddo
else
   !Simply compute y displacement for use in A equation:
  do ix=1,nx
    do iy=1,ny
      yy=gly*(y0(iy,ix)-yig(iy))
      ddy(iy,ix)=yy-elly*dble(int(yy*hlyi))
      qs(iy,ix)=qs(iy,ix)+beta*ddy(iy,ix)
    enddo
  enddo
endif

 !Spread y displacement diffusively and add to A:
call spread(ddy)
 !Add onto A:
do ix=1,nx
  do iy=1,ny
    aa(iy,ix)=aa(iy,ix)+b0*ddy(iy,ix)
  enddo
enddo

 !Integrate also qd, adiabatically; the source term is added just below:
call interpol(qdpre,qd,x0,y0)

 !Add PV source term to qd:
do ix=1,nx
  do iy=1,ny
     !Obtain old source (at time t) at the departure grid point using
     !bi-linear interpolation of sqpre:
    ix0=1+int(x0(iy,ix))
    pxc=dble(ix0)-x0(iy,ix)
    px=one-pxc
    ix1=ixp(ix0)

    iy0=1+int(y0(iy,ix))
    pyc=dble(iy0)-y0(iy,ix)
    py=one-pyc
    iy1=iyp(iy0)

    sqod=pyc*(pxc*sqpre(iy0,ix0)+px*sqpre(iy0,ix1)) &
      &  +py*(pxc*sqpre(iy1,ix0)+px*sqpre(iy1,ix1))
  
     !Integrate in time using trapezoidal rule (sq is the new source
     !at time t+dt) and add to qd:
    qd(iy,ix)=qd(iy,ix)+hfdt*(sqod+sq(iy,ix))
  enddo
enddo

!---------------------------------------------------------------------
 !Add any point vortices here to force the flow (they are converted
 !to gridded PV values and added to qd in spectral space):
if (stoch) then
   !Initialise array to uptake the vorticity of the vortices:
  do ix=1,nx
    do iy=1,ny
      sqstoch(iy,ix)=zero
    enddo
  enddo

  if (ivor .eq. 1) then
     !Add nvor vortices (as +/- arbitrarily separated pairs):
    nvor=2*nint(f12*(dnvor*t-totnvor))

     !The vortices are place randomly in the domain:
    dnx=dble(nx)
    dny=dble(ny)
    svor=-vorvor
    do k=1,nvor
      svor=-svor

       !Divide each vortex's circulation among the corners
       !of the local grid box (inverse bi-linear interpolation);
      xx=dnx*rand(0)
      ix0=1+int(xx)
      pxc=dble(ix0)-xx
      px=one-pxc
      ix1=ixp(ix0)

      yy=dny*rand(0)
      iy0=1+int(yy)
      pyc=dble(iy0)-yy
      py=one-pyc
      iy1=iyp(iy0)

      sqstoch(iy0,ix0)=sqstoch(iy0,ix0)+svor*pyc*pxc
      sqstoch(iy0,ix1)=sqstoch(iy0,ix1)+svor*pyc*px
      sqstoch(iy1,ix0)=sqstoch(iy1,ix0)+svor*py*pxc
      sqstoch(iy1,ix1)=sqstoch(iy1,ix1)+svor*py*px
    enddo

  else
     !Add nvor dipoles:
    nvor=nint(dnvor*t-totnvor)

     !The dipoles are place randomly in the domain with random
     !orientation:
    dnx=dble(nx)
    dny=dble(ny)
    do k=1,nvor

      theta=twopi*rand(0)
      cth=cos(theta)
      sth=sin(theta)

      xx=dnx*rand(0)
      ix0=1+int(xx)
      pxc=dble(ix0)-xx
      px=one-pxc
      ix1=ixp(ix0)

      yy=dny*rand(0)
      iy0=1+int(yy)
      pyc=dble(iy0)-yy
      py=one-pyc
      iy1=iyp(iy0)

      sqstoch(iy0,ix0)=sqstoch(iy0,ix0)-vorvor*(pyc*cth+pxc*sth)
      sqstoch(iy0,ix1)=sqstoch(iy0,ix1)+vorvor*(pyc*cth -px*sth)
      sqstoch(iy1,ix0)=sqstoch(iy1,ix0)+vorvor*(pxc*sth -py*cth)
      sqstoch(iy1,ix1)=sqstoch(iy1,ix1)+vorvor*( px*sth +py*cth)
    enddo

  endif

   !Total number of vortices/dipoles added so far:
  totnvor=totnvor+dble(nvor)

   !Finally add forcing to qd:
  do ix=1,nx
    do iy=1,ny
      qd(iy,ix)=qd(iy,ix)+sqstoch(iy,ix)
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine source(qqsrc)

! Gets the source terms for the PV evolution equation.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed array:
double precision:: qqsrc(ny,nx)

!-----------------------------------------------------------------
 !Compute curl of Lorentz force as PV source term:
call lorentz(aa,qqsrc)

 !Implement Ekman and/or thermal damping:
if (heating) then
  if (friction) then 
   !Use thermal and Ekman damping:
    do ix=1,nx
      do iy=1,ny
        qqsrc(iy,ix)=qqsrc(iy,ix)+therm*(pp(iy,ix)-ppeq(iy,ix)) &
               & -rekman*(qq(iy,ix)+kdsq*pp(iy,ix))
      enddo
    enddo
  else
    do ix=1,nx
      do iy=1,ny
        qqsrc(iy,ix)=qqsrc(iy,ix)+therm*(pp(iy,ix)-ppeq(iy,ix))
      enddo
    enddo
  endif
else
   !Only use Ekman damping
  do ix=1,nx
    do iy=1,ny
      qqsrc(iy,ix)=qqsrc(iy,ix)-rekman*(qq(iy,ix)+kdsq*pp(iy,ix))
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine inversion(ppflag)

! Inverts Laplace's operator on PV anomaly (PV - beta*y) to obtain 
! the streamfunction (pp) and the velocity (uu,vv) = (-dpp/dy,dpp/dx).
! pp is returned only if ppflag is true.

implicit double precision(a-h,o-z)
implicit integer(i-n)

logical:: ppflag

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qc)

 !Combine fields to update qq with full field:
call combine(qq,qc,qs,qd)

 !Invert PV to obtain velocity field: 
call main_invert(qq,uu,vv,pp,ppflag)      

return
end subroutine

!=======================================================================

subroutine adapt(igsave,icsave)

! Adapts the time step dt to ensure that it is less than or equal to
! the minimum of dtfac/max(|zeta|_max) and C*dx/max(|u|_max,|b|_max)
! where dx is the grid spacing, u is the vector velocity field, and
! b is the vector magnetic field.  C = cfl_max is specified below.

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision,parameter:: dtfac=pi/40.d0, cflmax=1.d0

!----------------------------------------------------------------------
 !Compute the gradients of A to get B = (B_0 + A_y, -A_x)
call magfield(aa,aax,aay)

 !Compute accurate advection time step:
umax=small
zzrms=zero
zzmax=small
do ix=1,nx
  do iy=1,ny
    usq=uu(iy,ix)**2+vv(iy,ix)**2
    bsq=aax(iy,ix)**2+(b0+aay(iy,ix))**2
    umax=max(umax,usq,bsq)
    zz=qq(iy,ix)+kdsq*pp(iy,ix)
    zzrms=zzrms+zz**2
    zzmax=max(zzmax,abs(zz))
     !Here, qq = gridded total PV anomaly (q - beta*y).
  enddo
enddo
umax=sqrt(umax)
zzrms=sqrt(zzrms*dsumi)
dtacc=min(glx*cflmax/umax,dtfac/max(zzmax,srwfm))
 !The restriction on the maximum Rossby wave frequency (srwfm)
 !ensures that the fastest Rossby wave frequency is resolved.

!---------------------------------------------------------------------
 !Choose a new time step: 
if (dt .gt. zero) then 
  dt=min(dtacc,dtmax)
  if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
else
   !Limit max timestep to a data save time
  dtmax=min(tgsave,tcsave)
  dt=min(dtacc,dtmax)
endif 
hfdt=dt/two

 !Increment the integral of max|zz|:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
 !Record various diagnostics to monitor.asc:
write(17,'(1x,f12.5,4(1x,f12.6),1x,f5.3)') & 
     & t,f12*zzrms**2,zzrms,zzmax,umax,twist

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
tgrid=tgsave*dble(itime)+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  igsave=1
   !Compute kinetic and magnetic energy:
  call binorm(pp,qq,eneupre)
  eneupre=-f12*eneupre
  call magfield(aa,aax,aay)
  call l2norm(aax,aaxl2)
  call l2norm(aay,aayl2)
  enebpre=f12*(aaxl2+aayl2)
else
   !Do not save data:
  igsave=0
endif

 !Set flag to save contour data every tcsave time units:
itime=int((t+dt)/tcsave)
tcont=tcsave*dble(itime)
if (t .le. tcont .and. t+dt .gt. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  icsave=sign(1,2*itime-1)
   !This construction avoids saving the contours at t = 0.
else
   !Do not save data:
  icsave=0
endif

return
end subroutine

!=======================================================================
      
subroutine savegrid

! Saves PV, energy and enstrophy spectrum at the desired save time

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: jj(ny,nx),jjl2
double precision:: qspec(0:max(nx,ny)),jspec(0:max(nx,ny))
real:: qqr4(ny,nx),tr4

 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Compute kinetic and magnetic energy:
call binorm(pp,qq,eneupost)
eneupost=-f12*eneupost
call magfield(aa,aax,aay)
call l2norm(aax,aaxl2)
call l2norm(aay,aayl2)
enebpost=f12*(aaxl2+aayl2)

 !Compute time interpolated energies:
eneu=pt*eneupre+ptc*eneupost
eneb=pt*enebpre+ptc*enebpost

 !Write energy to ene.asc:
write(15,'(f7.2,3(1x,f16.9))') tgrid,eneu,eneb,eneu+eneb

 !Interpolate PV anomaly and A at save time:
do ix=1,nx
  do iy=1,ny
    ula(iy,ix)=pt*qspre(iy,ix)+ptc*qq(iy,ix)
    vla(iy,ix)=pt*aapre(iy,ix)+ptc*aa(iy,ix)
  enddo
enddo

!Write full PV to qq.r4 and current density to jj.r4:
tr4=real(tgrid)
do ix=1,nx
  do iy=1,ny
    qqr4(iy,ix)=real(ula(iy,ix)+bety(iy))
  enddo
enddo
write(31,rec=igrids) tr4,qqr4

call current(vla,jj)
do ix=1,nx
  do iy=1,ny
    qqr4(iy,ix)=real(jj(iy,ix))
  enddo
enddo
write(32,rec=igrids) tr4,qqr4

 !Compute domain integral of (q-beta*y)^2 & j^2 and write to norms.asc:
call l2norm(ula,qql2)
call l2norm(jj,jjl2)
write(16,'(f7.2,2(1x,f16.9))') tgrid,qql2,jjl2

write(*,'(a,f7.2,4(a,f13.6))') ' t = ',tgrid, &
   & ' <q^2> = ',qql2,' <j^2> = ',jjl2,' E_u = ',eneu,' E_b = ',eneb

 !Compute 1d PV and current density spectra:
call spec1d(ula,qspec,0)
call spec1d(vla,jspec,1)
sumqspec=zero
sumjspec=zero
do k=1,kmax
  sumqspec=sumqspec+qspec(k)
  sumjspec=sumjspec+jspec(k)
   !Normalise to take into account uneven sampling of wavenumbers 
   !in each shell [k-1/2,k+1/2]:
  qspec(k)=spmf(k)*qspec(k)
  jspec(k)=spmf(k)*jspec(k)
enddo
sumqspec=8.d0*sumqspec*dsumi
sumjspec=8.d0*sumjspec*dsumi
 !Write out spectrum to file:
write(51,'(f7.2,4(1x,f16.9),1x,i5)') tgrid,sumqspec,qql2,sumjspec,jjl2,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(51,'(3(1x,f12.8))') alk(k),log10(qspec(k)),log10(jspec(k))
enddo
 !Note: alk(k) = log_10(k)

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves PV contours and residual PV for post-processing and imaging

implicit double precision(a-h,o-z)
implicit integer(i-n)

real:: qdr4(ny,nx),tr4
character(len=3):: pind

write(*,'(a,f12.5)') ' Saving contours at t = ',t
irec=nint(t/tcsave)
write(pind(1:3),'(i3.3)') irec

tr4=real(t)
 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
do ix=1,nx
  do iy=1,ny
    qdr4(iy,ix)=real(qq(iy,ix)-qc(iy,ix))
  enddo
enddo
write(83,rec=irec) tr4,qdr4

 !Save PV contours if any exist:
if (nq .gt. 0) then
  open(81,file='cont/qqindex'//pind,form='unformatted', &
      & access='direct',status='replace',recl=12*nq)
  write(81,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq)
  close(81)

  open(82,file='cont/qqnodes'//pind,form='unformatted', &
      & access='direct',status='replace',recl=16*nptq)
  write(82,rec=1) xq(1:nptq),yq(1:nptq)
  close(82)
endif

return
end subroutine

!=======================================================================

 !Main end module
end module
