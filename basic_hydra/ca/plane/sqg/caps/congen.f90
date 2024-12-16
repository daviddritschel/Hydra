module congen
! Converts contours (xq,yq) to gridded values on an ultra-fine grid 
! of dimensions mgu*nx x mgu*ny (in a closed rectangular domain), 
! optionally adds a residual field (interpolated to the ultra-fine 
! grid), and then creates new contours.

! Open contours originating and terminating in a boundary added by
! D G Dritschel on 18 June 2012 @ Moscow

use common

implicit none

double precision:: qa(nyu+1,nxu)
double precision:: xa(npm),ya(npm)
integer:: inda(nm),npa(nm),i1a(nm),i2a(nm)
integer:: na,npta

contains

!=====================================================================

subroutine recontour(qq)
! Main routine for recontouring (from D & Ambaum, 1996, QJRMS)
! qq: a gridded field added to that due to contours

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed array:
double precision:: qq(ny,nx)

 !-----------------------------------------------------------------
 !Counters for total number of nodes and contours:
npta=0
na=0

 !Obtain fine grid field qa:
if (nq .gt. 0) then
   !Convert contours to gridded values (in the array qa):
  call con2ugrid

   !Bi-linear interpolate the residual qq to the fine grid and add to qa:
  do ix=1,nxu
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    do iy=1,nyu
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq(iy0,ix0) &
                       & +w10(iyf,ixf)*qq(iy1,ix0) &
                       & +w01(iyf,ixf)*qq(iy0,ix1) &
                       & +w11(iyf,ixf)*qq(iy1,ix1)

    enddo
  enddo

else

   !Check if field requires contouring by computing l1 norm of qq:
  qql1=garea*sum(abs(qq))
  if (qql1 .lt. small) return

   !No contours: interpolate qq (which here contains the full field)
   !to the fine grid as qa:
  do ix=1,nxu
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    do iy=1,nyu
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qa(iy,ix)=w00(iyf,ixf)*qq(iy0,ix0)+w10(iyf,ixf)*qq(iy1,ix0) &
             & +w01(iyf,ixf)*qq(iy0,ix1)+w11(iyf,ixf)*qq(iy1,ix1)
 
    enddo
  enddo

endif

 !Generate new contours (xa,ya) from qa array:
call ugrid2con

 !Copy arrays back to those in the argument of the subroutine:
do i=1,npta
  xq(i)=xa(i)
  yq(i)=ya(i)
enddo

do j=1,na
  i1q(j)=i1a(j)
  i2q(j)=i2a(j)
  npq(j)=npa(j)
  indq(j)=inda(j)
enddo

nq=na
nptq=npta

return
end subroutine

!==========================================================================

subroutine ugrid2con
! Generates contours (xa,ya) from the gridded field qa for the levels
! +/-qjump/2, +/-3*qjump/2, ....

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local parameters and arrays:
integer,parameter:: ncrm=3*nplm/4
 !ncrm:  max number of contour crossings of a single contour level
 !nplm:  max number of nodes in any contour level
 
integer(kind=dbleint),parameter:: nxny=int(nxu,kind=dbleint)*int(nyu,kind=dbleint)
integer(kind=dbleint):: k,kob,kib(ncrm)
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qdx(nxu+1),qdy(nyu+1)
double precision:: xd(nprm),yd(nprm)
integer:: isx(nxu+1),isy(nyu+1)
integer:: icre(nm)
integer:: icrtab(nxny,2)
integer(kind=halfint):: noctab(nxny)
logical:: free(ncrm),keep

 !initialise constants and arrays:
qjumpi=one/qjump
qoff=qjump*dble(nlevm)
 !qoff: should be a large integer multiple of the contour interval, qjump.  
 !The multiple should exceed the maximum expected number of contour levels.

 !--------------------------------------------------------
 !First get the beginning and ending contour levels:
qamax=qa(1,1)
qamin=qa(1,1)
do ix=1,nxu
  do iy=1,nyu
    qamax=max(qamax,qa(iy,ix))
    qamin=min(qamin,qa(iy,ix))
  enddo
enddo

levbeg=int((qoff+qamin)*qjumpi+f12)+1
levend=int((qoff+qamax)*qjumpi+f12)

if (levbeg .le. levend) then
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Loop over contour levels and process:
do lev=levbeg,levend
 !Integer index giving contour level:
levt=lev-nlevm+(lev-1)/nlevm-1

 !Counter for total number of grid line crossings:
ncr=0

 !Contour level being sought:
qtmp=(dble(lev)-f12)*qjump-qoff

 !Below, kib = grid box into which the contour (containing ncr) is going
 !       kob =   "   "  out of "    "     "         "       "    " coming
 !      [kob -> ncr -> kib:  ncr lies at the boundary between kob & kib]

 !   *** grid boxes are numbered 1 (lower left) to nxu*nyu (upper right) ***

 !Initialise number of crossings per box:
do k=1,nxny
  noctab(k)=0
enddo

 !-----------------------------------------------------------
 !Find x grid line crossings first:
do ix=1,nxu
  xgt=xgu(ix)

  do iy=1,nyu
    qdy(iy)=qa(iy,ix)-qtmp
    isy(iy)=sign(one,qdy(iy))
  enddo
  qdy(nyu+1)=qdy(1)
  isy(nyu+1)=isy(1)

  do iy=1,nyu
    if (isy(iy) .ne. isy(iy+1)) then
      ncr=ncr+1
      inc=(1-isy(iy))/2
      kaa=(iy-1)*nxu+1
      kib(ncr)=kaa+ibx(ix+inc)
      kob=kaa+ibx(ix+1-inc)
      noctab(kob)=noctab(kob)+1
      icrtab(kob,noctab(kob))=ncr
      xcr(ncr)=xgt
      yy=ygu(iy)-glyu*qdy(iy)/(qdy(iy+1)-qdy(iy))
      ycr(ncr)=oms*(yy-elly*dble(int(yy*hlyi)))
    endif
  enddo

enddo

 !----------------------------------------------------------
 !Find y grid line crossings next:
do iy=1,nyu
  ygt=ygu(iy)

  do ix=1,nxu
    qdx(ix)=qa(iy,ix)-qtmp
    isx(ix)=sign(one,qdx(ix))
  enddo
  qdx(nxu+1)=qdx(1)
  isx(nxu+1)=isx(1)

  do ix=1,nxu
    if (isx(ix) .ne. isx(ix+1)) then
      ncr=ncr+1
      inc=(1-isx(ix))/2
      kib(ncr)=nxu*iby(iy+1-inc)+ix
      kob=nxu*iby(iy+inc)+ix
      noctab(kob)=noctab(kob)+1
      icrtab(kob,noctab(kob))=ncr
      ycr(ncr)=ygt
      xx=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
      xcr(ncr)=oms*(xx-ellx*dble(int(xx*hlxi)))
    endif
  enddo

enddo

 !----------------------------------------------------------------
 !Now re-build contours:
do icr=1,ncr
  free(icr)=.true.
enddo

do icr=1,ncr
  if (free(icr)) then
     !A new contour (indexed na) starts here:
    na=na+1
    inda(na)=levt
    ibeg=npta+1
    i1a(na)=ibeg

     !First point on the contour:
    npd=1
    xd(1)=xcr(icr)
    yd(1)=ycr(icr)

     !Find remaining points on the contour:
    k=kib(icr)
     !k is the box the contour is entering
    noc=noctab(k)
     !Use last crossing (noc) in this box (k) as the next node:
    icrn=icrtab(k,noc)
     !icrn gives the next point after icr (icrn is leaving box k)
    do while (icrn .ne. icr)
      noctab(k)=noc-1
       !noctab is usually zero now except for boxes with a
       !maximum possible 2 crossings
      npd=npd+1
      xd(npd)=xcr(icrn)
      yd(npd)=ycr(icrn)
      free(icrn)=.false.
      k=kib(icrn)
      noc=noctab(k)
      icrn=icrtab(k,noc)
    enddo

     !Re-distribute nodes on this contour 3 times to reduce complexity:
    keep=.false.
    do
      call renode(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
       !Delete contour if deemed too small (see renode_closed):
      if (npa(na) .eq. 0) exit
      call renode(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
       !Delete contour if deemed too small (see renode_closed):
      if (npd .eq. 0) exit
      call renode(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
       !Delete contour if deemed too small (see renode_closed):
      if (npa(na) .eq. 0) exit
       !Contour is big enough to keep:
      keep=.true.
      exit
    enddo

    if (keep) then 
      npta=npta+npa(na)
      iend=ibeg+npa(na)-1
      i2a(na)=iend
      do i=ibeg,iend-1
        nextq(i)=i+1
      enddo
      nextq(iend)=ibeg
    else
      na=na-1
    endif

    free(icr)=.false.
  endif
enddo

enddo
 !End of loop over contour levels
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
endif

return
end subroutine

!=======================================================================

subroutine con2ugrid
! Contour -> grid conversion.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local parameters and arrays:
double precision:: qjx(nxu)
double precision:: dx(nptq),dy(nptq)
integer:: ixc(nptq),nxc(nptq)
logical:: crossx(nptq)

!----------------------------------------------------------------
 !Initialise interior x grid line crossing information and fill the
 !q jump array along lower boundary:
do i=1,nptq
  ixc(i)=int(glxui*(xq(i)-xmin))
enddo

do ix=1,nxu
  qjx(ix)=zero
enddo

do i=1,nptq
  ia=nextq(i)
  xx=xq(ia)-xq(i)
  dx(i)=xx-ellx*dble(int(xx*hlxi))
  yy=yq(ia)-yq(i)
  dy(i)=yy-elly*dble(int(yy*hlyi))
  ixdif=ixc(ia)-ixc(i)
  nxc(i)=ixdif-nxu*((2*ixdif)/nxu)
  crossx(i)=(nxc(i) .ne. 0)
  if (abs(yy) .gt. hly) then
     !The contour segment (i,ia) crosses y = ymin; find x location:
    cc=sign(one,dy(i))
    p=-(yq(i)-hly*cc)/(dy(i)+small)
    ix=1+int(glxui*(mod(ellx+hlx+xq(i)+p*dx(i),ellx)))
    qjx(ix)=qjx(ix)-qjump*cc
     !Note: qjx gives the jump going from ix to ix+1
  endif
enddo

 !Sum q jumps to obtain the gridded q along lower boundary:
qa(1,1)=zero
 !Corner value cannot be determined a priori; qavg is used for this below
do ix=1,nxu-1
  qa(1,ix+1)=qa(1,ix)+qjx(ix)
enddo

!----------------------------------------------------------------
 !Initialise interior q jump array:
do ix=1,nxu
  do iy=2,nyu+1
    qa(iy,ix)=zero
  enddo
enddo

 !Determine x grid line crossings and accumulate q jumps:
do i=1,nptq
  if (crossx(i)) then
    jump=sign(1,nxc(i))
    ixbeg=ixc(i)+(1+jump)/2+nxu
    sqjump=qjump*jump
    ncr=0
    do while (ncr .ne. nxc(i)) 
      ix=1+mod(ixbeg+ncr,nxu)
      xx=xgu(ix)-xq(i)
      px0=(xx-ellx*dble(int(xx*hlxi)))/dx(i)
       !The contour crossed the fine grid line ix at the point
       !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
      yy=yq(i)+px0*dy(i)
      iy=2+int(glyui*(yy-elly*dble(int(yy*hlyi))-ymin))
       !Increment q jump between the grid lines iy-1 & iy:
      qa(iy,ix)=qa(iy,ix)+sqjump
       !Go on to consider next x grid line (if there is one):
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get q values by sweeping through y:
do ix=1,nxu
  do iy=2,nyu
    qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
  enddo
enddo

 !Remove average:
qavg=zero
do ix=1,nxu
  do iy=1,nyu
    qavg=qavg+qa(iy,ix)
  enddo
enddo
qavg=qavg/dble(nxu*nyu)

do ix=1,nxu
  do iy=1,nyu
    qa(iy,ix)=qa(iy,ix)-qavg
  enddo
enddo

return
end subroutine

!==========================================================================

 !Main end module
end module
