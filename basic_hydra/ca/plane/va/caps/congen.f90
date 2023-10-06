module congen

! Module containing subroutines for rebuilding PV contours given either
! the PV field and no contours, or the PV residual field together with
! contours.  This creates new contours in either case.

use common

implicit none

 !Array for storing the PV field interpolated to the ultra-fine grid: 
double precision:: qa(ngu+1,ngu)

 !Temporary arrays for contour storage:
double precision:: xa(npm),ya(npm)
integer:: inda(nm),npa(nm),i1a(nm),i2a(nm)
integer:: na,npta
 !Note: contours are built level by level to keep storage to a minimum.


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=====================================================================

subroutine recontour(qg)
! Main routine for recontouring (from Dritschel & Ambaum, 1996, QJRMS).
! qg: a gridded field added to that associated with any contours

implicit none

 !Passed array:
double precision:: qg(ng,ng)

 !Local variables:
integer:: ix,ixf,ix0,ix1
integer:: iy,iyf,iy0,iy1
integer:: i,j

 !-----------------------------------------------------------------
 !Counters for total number of nodes and contours:
npta=0
na=0

 !Obtain ultra-fine grid field qa:
if (nq .gt. 0) then
   !There are contours; convert them to gridded values (in the array qa):
  call con2ugrid

   !Bi-linear interpolate the residual qg to the fine grid and add to qa:
  do ix=1,ngu
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    do iy=1,ngu
      iyf=ixfw(iy)
      iy0=ix0w(iy)
      iy1=ix1w(iy)

      qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qg(iy0,ix0)+w10(iyf,ixf)*qg(iy1,ix0) &
                         +w01(iyf,ixf)*qg(iy0,ix1)+w11(iyf,ixf)*qg(iy1,ix1)

    enddo
  enddo

else
   !There are no contours (the usual situation at t = 0):

   !Check if field requires contouring by computing l1 norm of qg:
  if (garea*sum(abs(qg)) .lt. small) return

   !Interpolate qg (containing the full field) to the fine grid as qa:
  do ix=1,ngu
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    do iy=1,ngu
      iyf=ixfw(iy)
      iy0=ix0w(iy)
      iy1=ix1w(iy)

      qa(iy,ix)=w00(iyf,ixf)*qg(iy0,ix0)+w10(iyf,ixf)*qg(iy1,ix0) &
               +w01(iyf,ixf)*qg(iy0,ix1)+w11(iyf,ixf)*qg(iy1,ix1)
 
    enddo
  enddo

endif

 !Generate new contours (xa,ya) from the gridded field in qa:
call ugrid2con

 !Copy contour arrays and indices back to those in the main code:
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

implicit none

 !Local parameters and arrays:
integer,parameter:: ncrm=3*nplm/4
 !ncrm:  maximum number of contour crossings of a single contour level
 !nplm:  maximum number of nodes in any contour level
 
integer(kind=dbleint),parameter:: ngsq=int(ngu,kind=dbleint)*int(ngu,kind=dbleint)
integer(kind=dbleint):: k,kob,kaa,kib(ncrm)
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qdx(ngu+1),qdy(ngu+1)
double precision:: xd(nprm),yd(nprm)
double precision:: qoff,qjumpi,sqjump,qtmp,xgt,ygt,xx,yy
integer:: isx(ngu+1),isy(ngu+1),icre(nm),icrtab(ngsq,2)
integer:: levbeg,levend,lev,levt,noc,icrn
integer:: ncr,i,ix,iy,icr,inc,ibeg,iend,npd
integer(kind=halfint):: noctab(ngsq)
logical:: free(ncrm),keep

 !----------------------------------------
 !Initialise various constants used below:
qjumpi=one/qjump
qoff=qjump*dble(nlevm)
 !qoff: should be a large integer multiple of the contour interval, qjump.  
 !The multiple should exceed the maximum expected number of contour levels.

 !--------------------------------------------------
 !First get the beginning and ending contour levels:
levbeg=int((qoff+minval(qa(1:ngu,:)))*qjumpi+f12)+1
levend=int((qoff+maxval(qa(1:ngu,:)))*qjumpi+f12)

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

   !*** grid boxes are numbered 1 (lower left) to ngu*ngu (upper right) ***

   !Initialise number of crossings per box:
  noctab=0

   !-----------------------------------------------------------
   !Find x grid line crossings first:
  do ix=1,ngu
    xgt=xgu(ix)

    qdy(1:ngu)=qa(1:ngu,ix)-qtmp
    isy(1:ngu)=sign(one,qdy(1:ngu))
    qdy(ngu+1)=qdy(1)
    isy(ngu+1)=isy(1)

    do iy=1,ngu
      if (isy(iy) .ne. isy(iy+1)) then
        ncr=ncr+1
        inc=(1-isy(iy))/2
        kaa=(iy-1)*ngu+1
        kib(ncr)=kaa+ibx(ix+inc)
        kob=kaa+ibx(ix+1-inc)
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        xcr(ncr)=xgt
        yy=xgu(iy)-glu*qdy(iy)/(qdy(iy+1)-qdy(iy))
        ycr(ncr)=oms*(yy-twopi*dble(int(yy*pinv)))
      endif
    enddo

  enddo

   !----------------------------------------------------------
   !Find y grid line crossings next:
  do iy=1,ngu
    ygt=xgu(iy)

    qdx(1:ngu)=qa(iy,1:ngu)-qtmp
    isx(1:ngu)=sign(one,qdx(1:ngu))
    qdx(ngu+1)=qdx(1)
    isx(ngu+1)=isx(1)

    do ix=1,ngu
      if (isx(ix) .ne. isx(ix+1)) then
        ncr=ncr+1
        inc=(1-isx(ix))/2
        kib(ncr)=ngu*ibx(iy+1-inc)+ix
        kob=ngu*ibx(iy+inc)+ix
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        ycr(ncr)=ygt
        xx=xgu(ix)-glu*qdx(ix)/(qdx(ix+1)-qdx(ix))
        xcr(ncr)=oms*(xx-twopi*dble(int(xx*pinv)))
      endif
    enddo

  enddo

   !----------------------------------------------------------------
   !Now re-build contours:
  free(1:ncr)=.true.

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
         !Delete contour if deemed too small (see renode in contours.f90):
        if (npa(na) .eq. 0) exit
        call renode(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
         !Delete contour if deemed too small:
        if (npd .eq. 0) exit
        call renode(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
         !Delete contour if deemed too small:
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

implicit none

 !Local parameters and arrays:
double precision:: qjx(ngu)
double precision:: dx(nptq),dy(nptq)
double precision:: xx,yy,cc,p,sqjump,px0,qavg
integer:: ixc(nptq),ngc(nptq)
integer:: i,ia,ix,ixdif,ixbeg,iy,jump,ncr
logical:: crossx(nptq)

!-------------------------------------------------------------------
 !Initialise interior x grid line crossing information and fill the
 !PV (q) jump array along the lower boundary (iy = 1 or y = -pi):
do i=1,nptq
  ixc(i)=int(glui*(xq(i)+pi))
enddo

qjx=zero

do i=1,nptq
  ia=nextq(i)
  xx=xq(ia)-xq(i)
  dx(i)=xx-twopi*dble(int(xx*pinv))
  yy=yq(ia)-yq(i)
  dy(i)=yy-twopi*dble(int(yy*pinv))
  ixdif=ixc(ia)-ixc(i)
  ngc(i)=ixdif-ngu*((2*ixdif)/ngu)
  crossx(i)=(ngc(i) .ne. 0)
  if (abs(yy) .gt. pi) then
     !The contour segment (i,ia) crosses y = -pi; find x location:
    cc=sign(one,dy(i))
    p=-(yq(i)-pi*cc)/(dy(i)+small)
    ix=1+int(glui*(mod(thrpi+xq(i)+p*dx(i),twopi)))
    qjx(ix)=qjx(ix)-qjump*cc
     !Note: qjx gives the jump in q going from ix to ix+1
  endif
enddo

 !Sum q jumps to obtain the gridded field qa along lower boundary:
qa(1,1)=zero
 !Corner value cannot be determined a priori; qavg is used for this below
do ix=1,ngu-1
  qa(1,ix+1)=qa(1,ix)+qjx(ix)
enddo

!----------------------------------------------------------------
 !Initialise interior q jump array:
qa(2:ngu+1,:)=zero

 !Determine x grid line crossings and accumulate q jumps:
do i=1,nptq
  if (crossx(i)) then
    jump=sign(1,ngc(i))
    ixbeg=ixc(i)+(1+jump)/2+ngu
    sqjump=qjump*dble(jump)
    ncr=0
    do while (ncr .ne. ngc(i)) 
      ix=1+mod(ixbeg+ncr,ngu)
      xx=xgu(ix)-xq(i)
      px0=(xx-twopi*dble(int(xx*pinv)))/dx(i)
       !The contour crossed the fine grid line ix at the point
       !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
      yy=yq(i)+px0*dy(i)
      iy=2+int(glui*(yy-twopi*dble(int(yy*pinv))+pi))
       !Increment q jump between the grid lines iy-1 & iy:
      qa(iy,ix)=qa(iy,ix)+sqjump
       !Go on to consider next x grid line (if there is one):
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get q values by sweeping through y:
do ix=1,ngu
  do iy=2,ngu
    qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
  enddo
enddo

 !Remove domain average:
qavg=sum(qa(1:ngu,:))/dble(ngu*ngu)
qa(1:ngu,:)=qa(1:ngu,:)-qavg

return
end subroutine

!==========================================================================

 !Main end module
end module
