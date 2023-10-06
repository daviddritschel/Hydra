module congen

!---------------------------------------------------------------------------
!         Converts PV contours to gridded values on an
!         ultra-fine grid of dimensions mgu*nt x mgu*ng (periodic 
!         in x but free slip walls in y), adds the residual q 
!         (interpolated to the ultra-fine grid), then creates new contours.  

!         This routine processes one contour level at a time.

!         This avoids storing large arrays in memory.  After
!         running this programme.

!         Adapted from ~dgd/cs/spe/sources/mcongen.F on 26/2/13
!         by Stuart King & DG Dritschel @ St Andrews

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

use common

implicit none

 !Grid -> Contour arrays:
double precision:: qa(0:ngu+1,ntu+1)

contains
 
!==========================================================================

subroutine recontour

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local parameters and variables:
double precision:: qq0(0:ng,nt)
double precision:: xtd(nt),utd(nt)

!==========================================================
 !Loop over great circles in latitude (only half the longitudes):
do ix=1,ng
  ic=ix+ng

   !Source vector:
  utd(1)=f23*(qd(1,ix)+qd(1,ic))
  do j=2,ng
    utd(j)=f23*(qd(j,ix)+qd(j-1,ix))
  enddo
  utd(ngp1)=f23*(qd(ng,ic)+qd(ng,ix))
  do j=ngp2,nt
    utd(j)=f23*(qd(ntp2-j,ic)+qd(ntp1-j,ic))
  enddo

   !Interpolate qd by 4th-order method (periodic):
  xtd(1)=utd(1)*htd(1)
  do j=2,nt
    xtd(j)=(utd(j)-f16*xtd(j-1))*htd(j)
  enddo
  do j=ntm2,1,-1
    xtd(j)=etd(j)*xtd(j+1)+xtd(j)
  enddo
  xtd(nt)=(etd(nt)*xtd(1)+xtd(nt))*xndeno
  xend=xtd(nt)

  do j=1,ntm1
    xtd(j)=ptd(j)*xend+xtd(j)
  enddo

   !Copy back into full grid array (qq0):
  do j=0,ng
    qq0(j,ix)=xtd(j+1)
  enddo
  qq0(0,ic)=xtd(1)
  do j=1,ng
    qq0(j,ic)=xtd(ntp1-j)
  enddo

enddo
 !Ends loops over great circles.  Interpolation complete.

 !Obtain unique polar values of qd for use below:
qdsp=zero
qdnp=zero
do ix=1,nt
  qdsp=qdsp+qq0(0 ,ix)
  qdnp=qdnp+qq0(ng,ix)
enddo
qdsp=qdsp/dble(nt)
qdnp=qdnp/dble(nt)
do ix=1,nt
  qq0(0 ,ix)=qdsp
  qq0(ng,ix)=qdnp
enddo

!------------------------------------------------------------
 !Obtain gridded PV from contours (if present):
if (n .gt. 0) then
   !Convert contours to gridded values:
  call con2ugrid

   !Bi-linear interpolate qq0 to the fine grid and add to qa:
  do ix=1,ntu
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)
  
    qa(0,ix)=qa(0,ix)+qdsp
    do iy=1,ngu-1
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq0(iy0,ix0)+w10(iyf,ixf)*qq0(iy1,ix0) &
                       & +w01(iyf,ixf)*qq0(iy0,ix1)+w11(iyf,ixf)*qq0(iy1,ix1)
    enddo
    qa(ngu,ix)=qa(ngu,ix)+qdnp
  enddo
   !qdsp & qdnp are the polar qd values (necessarily uniform).

else
   !No contours: interpolate qq0 to the fine grid as qa:
  do ix=1,ntu
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    qa(0,ix)=qdsp
    do iy=1,ngu-1
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)
      qa(iy,ix)=w00(iyf,ixf)*qq0(iy0,ix0)+w10(iyf,ixf)*qq0(iy1,ix0) &
               +w01(iyf,ixf)*qq0(iy0,ix1)+w11(iyf,ixf)*qq0(iy1,ix1)
    enddo
    qa(ngu,ix)=qdnp
  enddo
   !qdsp & qdnp are the polar q values (necessarily uniform).

endif

 !Remove global average:
avqa=zero
do ix=1,ntu
  do iy=1,ngu-1
    avqa=avqa+qa(iy,ix)*rdtu(iy)
  enddo
enddo
 !rdtu = rho/tau on the fine grid (see contours.f90)
avqa=avqa*dsumui

do ix=1,ntu
  do iy=0,ngu
    qa(iy,ix)=qa(iy,ix)-avqa
  enddo
enddo

 !Add a periodic column at ix = ntu+1:
ix=ntu+1
do iy=0,ngu
  qa(iy,ix)=qa(iy,1)
enddo

 !Counters for total number of nodes and contours:
npt=0
n=0

 !Generate new contours:
call ugrid2con

return
end subroutine

!=======================================================================

subroutine con2ugrid
! Converts PV contours (x,y,z) to gridded values (qa).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local arrays:
double precision:: cx(npt),cy(npt),cz(npt)
double precision:: sq(npt)
integer:: ntc(npt),ilm1(npt)

!----------------------------------------------------------------
 !Initialise crossing information:
do k=1,npt
  ilm1(k)=int(dlui*(pi+atan2(y(k),x(k))))
enddo

do k=1,npt
  ka=next(k)
  cx(k)=z(k)*y(ka)-y(k)*z(ka)
  cy(k)=x(k)*z(ka)-z(k)*x(ka)
  cz(k)=x(k)*y(ka)-y(k)*x(ka)
  ntc(k)=ilm1(ka)-ilm1(k)
enddo

do k=1,npt
  sig=sign(one,cz(k))
  sq(k)=dq*sig
  ntc(k)=ntc(k)-ntu*((2*ntc(k))/ntu)
  if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
  if (abs(cz(k)) .gt. zero) then
    cx(k)=cx(k)/cz(k)
    cy(k)=cy(k)/cz(k)
  endif
enddo

!----------------------------------------------------------------------
 !Initialise PV jump array:
do i=1,ntu
  do j=0,ngu+1
    qa(j,i)=zero
  enddo
enddo

 !Determine crossing indices:
do k=1,npt
  if (ntc(k) .ne. 0) then
    jump=sign(1,ntc(k))
    ioff=ntu+ilm1(k)+(1+jump)/2
    ncr=0
    do while (ncr .ne. ntc(k))
      i=1+mod(ioff+ncr,ntu)
      rlatc=dlui*(hpi+atan(cx(k)*clonu(i)+cy(k)*slonu(i)))
      j=int(rlatc)+1
      p=rlatc-dble(j-1)
      qa(j,i)=  qa(j,i)+(one-p)*sq(k)
      qa(j+1,i)=qa(j+1,i)+    p*sq(k)
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get PV values, at half latitudes, by sweeping through latitudes:
do i=1,ntu
  do j=2,ngu
    qa(j,i)=qa(j,i)+qa(j-1,i)
  enddo
enddo
 !Here, qa(j,i) stands for the PV at latitude j-1/2,
 !from j = 1, ..., ngu.

 !Determine unique polar values:
qasp=zero
qanp=zero
do i=1,ntu
  qasp=qasp+qa(1  ,i)
  qanp=qanp+qa(ngu,i)
enddo
qasp=qasp/dble(ntu)
qanp=qanp/dble(ntu)

 !Average half-grid PV to full grid:
do i=1,ntu
  qa(0,i)=qasp
  do j=1,ngu-1
    qa(j,i)=f12*(qa(j,i)+qa(j+1,i))
  enddo
  qa(ngu,i)=qanp
enddo

return
end subroutine

!=======================================================================
subroutine ugrid2con
! Generates new contours (xd,yd) from the gridded data (qa).

implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: nreno=2
! nreno: number of times renode is called to reduce point
!        density on contours.

integer,parameter:: ncrm=3*npm/2
! ncrm: max number of contour crossings of a single field
!       level on the finest grid

 !Local Grid -> Contour arrays:
integer(kind=dbleint),parameter:: ntng=int(ntu,kind=dbleint)*int(ngu,kind=dbleint)
integer(kind=dbleint):: kob,kib(ncrm)
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qdx(ntu+1),qdy(0:ngu)
double precision:: xd(nprm),yd(nprm),zd(nprm)
integer:: isx(ntu+1),isy(0:ngu)
integer:: icrtab(ntng,2)
integer(kind=halfint):: noctab(ntng)
logical:: free(ncrm),keep

!--------------------------------------------------------
 !First get the beginning and ending contour levels:
qamax=max(qa(0,1),qa(ngu,1))
qamin=min(qa(0,1),qa(ngu,1))
 !(qa is uniform at the extended edges iy = 0 and ngu)
do ix=1,ntu
  do iy=1,ngu-1
    qamax=max(qamax,qa(iy,ix))
    qamin=min(qamin,qa(iy,ix))
  enddo
enddo

levbeg=int((qoff+qamin)*dqi+f12)+1
levend=int((qoff+qamax)*dqi+f12)

 !Return if no levels to process:
if (levbeg .gt. levend) return 

 !Loop over contour levels and process:
do lev=levbeg,levend
   !Integer index giving contour level:                                                        
  indq=lev-nlevm+(lev-1)/nlevm-1
   !Counter for total number of grid line crossings:
  ncr=0

   !Contour level being sought:
  qtmp=qlev(lev)

   !Initialise number of crossings per box:
  do kob=1,ntng
    noctab(kob)=0
  enddo

   !Find x grid line crossings first:
  do ix=1,ntu
    xgt=xgu(ix)

    do iy=0,ngu
      qdy(iy)=qa(iy,ix)-qtmp
      isy(iy)=sign(one,qdy(iy))
    enddo

    do iy=0,ngu-1
      if (isy(iy) .ne. isy(iy+1)) then
        ncr=ncr+1
        inc=(1-isy(iy))/2
        kib(ncr)=iy+1+ibx(ix,inc)
        kob=iy+1+ibx(ix,1-inc)
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        xcr(ncr)=xgt
        ycr(ncr)=ygu(iy)-glyu*qdy(iy)/(qdy(iy+1)-qdy(iy))
      endif
    enddo

  enddo

!   Above, kib = grid box into which the contour (containing icr) is going
!          kob =   "   "  out of "    "     "         "       "    " coming
!     [kob -> icr -> kib:  icr lies at the boundary between kob & kib]

 !Find y grid line crossings next (no crossings can occur at iy=0,ngu):
  do iy=1,ngu-1
    ygt=ygu(iy)

    do ix=1,ntu+1
      qdx(ix)=qa(iy,ix)-qtmp
      isx(ix)=sign(one,qdx(ix))
    enddo

    do ix=1,ntu
      if (isx(ix) .ne. isx(ix+1)) then
        ncr=ncr+1
        inc=(1-isx(ix))/2
        kib(ncr)=ibx(ix,1)+iy+1-inc
        kob=ibx(ix,1)+iy+inc
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        ycr(ncr)=ygt
        xcr(ncr)=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
      endif
    enddo

  enddo

!------------------------------------------------------------------------
   !Now re-build contours - converting to spherical geometry:
  do icr=1,ncr
    free(icr)=.true.
  enddo

  do icr=1,ncr
    if (free(icr)) then
       !A new contour (indexed n) starts here:
      n=n+1
      ind(n)=indq
      ibeg=npt+1
      i1(n)=ibeg

       !First point on the contour:
      npd=1
      coslat=cos(ycr(icr))
      xd(1)=coslat*cos(xcr(icr))
      yd(1)=coslat*sin(xcr(icr))
      zd(1)=sin(ycr(icr))

       !Find remaining points on the contour:
      kob=kib(icr)
       !kib(icr) is the box the contour is entering
      noc=noctab(kob)
       !Use last crossing (noc) in this box (kob) as the next node:
      icrn=icrtab(kob,noc)
       !icrn gives the next point after icr (icrn is leaving box kob)
      do while (icrn .ne. icr)
        noctab(kob)=noc-1
         !noctab is usually zero now except for boxes with a
         !maximum possible 2 crossings
        npd=npd+1
        coslat=cos(ycr(icrn))
        xd(npd)=coslat*cos(xcr(icrn))
        yd(npd)=coslat*sin(xcr(icrn))
        zd(npd)=sin(ycr(icrn))
        free(icrn)=.false.
        kob=kib(icrn)
        noc=noctab(kob)
        icrn=icrtab(kob,noc)
      enddo

       !Re-distribute nodes on this contour 3 times to reduce complexity:
      keep=.false.
      do
        call renode(xd,yd,zd,npd,x(ibeg),y(ibeg),z(ibeg),np(n))
         !Delete contour if deemed too small (see renode):
        if (np(n) .eq. 0) exit
        call renode(x(ibeg),y(ibeg),z(ibeg),np(n),xd,yd,zd,npd)
         !Delete contour if deemed too small (see renode):
        if (npd .eq. 0) exit
        call renode(xd,yd,zd,npd,x(ibeg),y(ibeg),z(ibeg),np(n))
         !Delete contour if deemed too small (see renode):
        if (np(n) .eq. 0) exit
         !Contour is big enough to keep:
        keep=.true.
        exit
      enddo

      if (keep) then 
        npt=npt+np(n)
        iend=ibeg+np(n)-1
        i2(n)=iend
        do i=ibeg,iend-1
          next(i)=i+1
        enddo
        next(iend)=ibeg
      else
        n=n-1
      endif

      free(icr)=.false.
    endif
  enddo

enddo
!End of loop over contour levels
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

return
end subroutine

!=======================================================================

 !Main end module
end module
