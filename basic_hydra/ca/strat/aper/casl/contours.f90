module contours

! Module contains all subroutines related to contour advection, for the 
! casl suite of f90 codes.

use constants
use variables
use generic

implicit none

 !Vorticity contours:
double precision:: xz(npm),yz(npm),zjump,zavg
integer:: indz(nm),npz(nm),i1z(nm),i2z(nm)
integer:: nextz(npm),nptz,nz

 !Buoyancy contours:
double precision:: xb(npm),yb(npm),bjump,bavg
integer:: indb(nm),npb(nm),i1b(nm),i2b(nm)
integer:: nextb(npm),nptb,nb

 !Contour to grid conversion quantities (fine grid used in this module):
double precision,parameter:: dxxf =two*hlx/dble(nxf),dxxfi=dble(nxf)/(two*hlx)
double precision,parameter:: dyyf =two*hly/dble(nyf),dyyfi=dble(nyf)/(two*hly)
double precision:: xxf(0:nxf)

 !Fine grid needed for potential energy calculation (in this module):
double precision:: xgf(0:nxf),ygf(0:nyf)
double precision,parameter:: glxf =ellx/dble(nxf),glxfi=dble(nxf)/ellx
double precision,parameter:: glyf =elly/dble(nyf),glyfi=dble(nyf)/elly

 !Contour to grid conversion quantities (ultra-fine grid used in congen):
double precision,parameter:: dxxu =two*hlx/dble(nxu),dxxui=dble(nxu)/(two*hlx)
double precision,parameter:: dyyu =two*hly/dble(nyu),dyyui=dble(nyu)/(two*hly)
double precision:: xxu(0:nxu)

 !Ultra-fine grid needed for contour regeneration (module congen):
double precision:: xgu(0:nxu),ygu(0:nyu)
double precision,parameter:: glxu =ellx/dble(nxu),glxui=dble(nxu)/ellx
double precision,parameter:: glyu =elly/dble(nyu),glyui=dble(nyu)/elly

 !For computing db/dx from contours (in routine getzzsrc below):
double precision:: div(ndiv)
integer:: ixfp(0:nxf),iyfp(0:nyf)

 !Next grid points used in bilinear interpolation:
integer:: ixp(0:nx),iyp(0:ny)

 !Area weights used for interpolation in module congen:
double precision:: w00(mgu,mgu),w10(mgu,mgu),w01(mgu,mgu),w11(mgu,mgu)
integer:: ixfw(0:nxu),ix0w(0:nxu),ix1w(0:nxu)
integer:: iyfw(0:nyu),iy0w(0:nyu),iy1w(0:nyu)

 !Reference state potential energy and available potential energy:
double precision:: peref,ape
 
 !Fine-grid area weightings used in PE and vorticity tendency calculations:
double precision,parameter:: gareaf=glxf*glyf,dsumfi=one/dble(nxf*nyf)
double precision:: wdzdt

 !Surgery & node redistribution quantities:
double precision,parameter:: amu=0.2d0,ell=6.25d0*glx
double precision,parameter:: dm=amu**2*ell/four,dm2=dm**2,d4small=small*glx**4
double precision,parameter:: dmsq=four*dm2,dmi=two/dm
double precision,parameter:: elf=one/ell**2,densf=one/(amu*sqrt(ell))

contains 

!=======================================================================

subroutine init_contours

! Initialises all quantities needed for contour advection.

implicit double precision(a-h,o-z)
implicit integer(i-n)

!--------------------------------------------------------------------
 !Next grid points used in velocity interpolation (velint) and elsewhere; 
 !these are needed to avoid accessing array values at ix = nx+1 and iy = ny+1:
do ix=0,nxm1
  ixp(ix)=ix+1
enddo
ixp(nx)=nx

do iy=0,nym1
  iyp(iy)=iy+1
enddo
iyp(ny)=ny

 !Used for interpolation in routine getzzsrc below:
do ix=0,nxf-1
  ixfp(ix)=ix+1
enddo
ixfp(nxf)=nxf

do iy=0,nyf-1
  iyfp(iy)=iy+1
enddo
iyfp(nyf)=nyf

 !Fine dilated x-grid lines needed for contour-to-grid conversion:
do ix=0,nxf
  xxf(ix)=xbeg+dxxf*dble(ix)
enddo

 !Actual fine grid lines (offset by xmin, ymin) needed for 
 !potential energy calculation:
do ix=0,nxf
  xgf(ix)=glxf*dble(ix)-xcen
enddo

do iy=0,nyf
  ygf(iy)=glyf*dble(iy)-ycen
enddo

 !Divisions used for obtaining db/dx from contours:
fac=one/dble(ndiv)
do k=1,ndiv
  div(k)=fac*dble(k)
enddo
 !Weighting factor for transferring db/dx to gridded values:
wdzdt=bjump/gareaf

!=========================================================================
 !Initialise fixed quantities needed in contour regeneration (module congen)

 !Ultra-fine dilated x-grid lines needed for contour-to-grid conversion:
do ix=0,nxu
  xxu(ix)=xbeg+dxxu*dble(ix)
enddo

!---------------------------------------------------------------
 !Actual ultra-fine grid lines needed for contour regeneration:
do ix=0,nxu
  xgu(ix)=xmin+glxu*dble(ix)
enddo

do iy=0,nyu
  ygu(iy)=ymin+glyu*dble(iy)
enddo

!-------------------------------------------------------------
 !Initialise area weights for interpolation of a gridded field
 !onto the ultra-fine horizontal grid:
fac=one/dble(mgu)

do ixf=1,mgu
  pxf=fac*dble(ixf-1)
  pxc=one-pxf

  do iyf=1,mgu
    pyf=fac*dble(iyf-1)
    pyc=one-pyf

    w00(iyf,ixf)=pyc*pxc
    w10(iyf,ixf)=pyf*pxc
    w01(iyf,ixf)=pyc*pxf
    w11(iyf,ixf)=pyf*pxf
  enddo
enddo

do ix=0,nxu-1
  ixx=ix/mgu
  ix0w(ix)=ixx
  ix1w(ix)=1+ixx
  ixfw(ix)=1+ix-mgu*ixx
enddo
ix0w(nxu)=nx
ix1w(nxu)=nx
ixfw(nxu)=1

do iy=0,nyu-1
  iyy=iy/mgu
  iy0w(iy)=iyy
  iy1w(iy)=1+iyy
  iyfw(iy)=1+iy-mgu*iyy
enddo
iy0w(nyu)=ny
iy1w(nyu)=ny
iyfw(nyu)=1

return
end subroutine

!=======================================================================

subroutine getzzsrc(dzdt,t0)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local arrays:
double precision:: a(nptb),b(nptb),c(nptb),e(nptb)
double precision:: u(nptb),v(nptb),dx(nptb),dy(nptb)
double precision:: dzdt(0:ny,0:nx)
double precision:: dzdtf(0:nyf+1,0:nxf)

!--------------------------------------------------------------------
if (osci) then 
   !Domain is oscillating, obtain angle:
  arg=omegag*t0
  theta=theini+atilt*cos(arg)+btilt*sin(arg)
  ctheta=cos(theta)
  stheta=sin(theta)
endif

!--------------------------------------------------------------------
 !Initialise dzdtf to zero everywhere:
do ix=0,nxf
  do iy=0,nyf
    dzdtf(iy,ix)=zero
  enddo
enddo

 !Compute cubic interpolation coefficients:
do j=1,nb
  is=i1b(j)
  ie=i2b(j)

  if (nextb(ie) .eq. 0) then
     !Contour j is open; it starts and ends at an edge
    ie=ie-1
    do i=is,ie
      dx(i)=xb(i+1)-xb(i)
      a(i+1)=-dx(i)
      dy(i)=yb(i+1)-yb(i)
      c(i+1)=-dy(i)
      v(i)=dx(i)*dx(i)+dy(i)*dy(i)
      u(i+1)=v(i)
      e(i)=sqrt(v(i))
    enddo
    b(is)=zero
    do i=is+1,ie
      if (dx(i)*a(i)+dy(i)*c(i) .gt. zero) then
         !Set curvature to zero at corners:
        b(i)=zero
      else
        b(i)=(dx(i)*c(i)-a(i)*dy(i))/ &
         & sqrt((a(i)*v(i)-dx(i)*u(i))**2+(c(i)*v(i)-dy(i)*u(i))**2+small3)
      endif
    enddo
    b(ie+1)=zero
    do i=is,ie
      bsum=e(i)*(b(i+1)+b(i))
      bdif=e(i)*(b(i+1)-b(i))
      a(i)=f16*bdif-f12*bsum
      b(i)=f12*(bsum-bdif)
      c(i)=f13*bdif
    enddo
  else
     !Contour j is closed
    do i=is,ie
      ia=nextb(i)
      dx(i)=xb(ia)-xb(i)
      dy(i)=yb(ia)-yb(i)
      v(i)=dx(i)*dx(i)+dy(i)*dy(i)
      e(i)=sqrt(v(i))
    enddo
    do ib=is,ie
      i=nextb(ib)
      u(i)=v(ib)
      a(i)=-dx(ib)
      c(i)=-dy(ib)
    enddo
    do i=is,ie
      if (dx(i)*a(i)+dy(i)*c(i) .gt. zero) then
         !Set curvature to zero at corners:
        b(i)=zero
      else
        b(i)=(dx(i)*c(i)-a(i)*dy(i))/ &
         & sqrt((a(i)*v(i)-dx(i)*u(i))**2+(c(i)*v(i)-dy(i)*u(i))**2+small3)
      endif
    enddo
    do i=is,ie
      ia=nextb(i)
      bsum=e(i)*(b(ia)+b(i))
      bdif=e(i)*(b(ia)-b(i))
      a(i)=f16*bdif-f12*bsum
      b(i)=f12*(bsum-bdif)
      c(i)=f13*bdif
    enddo
  endif

  do i=is,ie
    x1=xb(i)
    y1=yb(i)
    do k=1,ndiv
      eta=div(k)*(a(i)+div(k)*(b(i)+div(k)*c(i)))
      x2=xb(i)+div(k)*dx(i)-eta*dy(i)
      y2=yb(i)+div(k)*dy(i)+eta*dx(i)
      avgsrc=wdzdt*((y1-y2)*ctheta+(x1-x2)*stheta)

      xx=glxfi*(f12*(x1+x2)-xmin)
      ix0=int(xx)
      ix1=ixfp(ix0)
      px=xx-dble(ix0)
      pxc=one-px

      yy=glyfi*(f12*(y1+y2)-ymin)
      iy0=int(yy)
      iy1=iyfp(iy0)
      py=yy-dble(iy0)
      pyc=one-py
       
      dzdtf(iy0,ix0)=dzdtf(iy0,ix0)+pyc*pxc*avgsrc
      dzdtf(iy0,ix1)=dzdtf(iy0,ix1)+pyc*px*avgsrc
      dzdtf(iy1,ix0)=dzdtf(iy1,ix0)+py*pxc*avgsrc
      dzdtf(iy1,ix1)=dzdtf(iy1,ix1)+py*px*avgsrc

      x1=x2
      y1=y2
    enddo
  enddo
enddo

 !Double the edge values (by symmetry of b):
do ix=0,nxf
  dzdtf(0,  ix)=two*dzdtf(0,  ix)
  dzdtf(nyf,ix)=two*dzdtf(nyf,ix)  
enddo
do iy=0,nyf
  dzdtf(iy,  0)=two*dzdtf(iy,  0)
  dzdtf(iy,nxf)=two*dzdtf(iy,nxf)
enddo

 !Average to inversion grid by repeated 1-2-1 averages in each direction:
call coarsen(dzdtf,dzdt)

return
end subroutine

!=======================================================================
      
subroutine velint(uu,vv,xq,yq,uq,vq,nptq)

! Bi-linearly interpolates current velocity (uu,vv) to the contour
! nodes (xq,yq) and stores the result in (uq,vq).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: uu(0:ny,0:nx),vv(0:ny,0:nx)
double precision:: xq(npm),yq(npm),uq(nptq),vq(nptq)

do i=1,nptq
  xx=glxi*(xq(i)-xmin)
  ix0=int(xx)
  ix1=ixp(ix0)
  px=xx-dble(ix0)
  pxc=one-px

  yy=glyi*(yq(i)-ymin)
  iy0=int(yy)
  iy1=iyp(iy0)
  py=yy-dble(iy0)
  pyc=one-py

  uq(i)=pyc*(pxc*uu(iy0,ix0)+px*uu(iy0,ix1)) &
      & +py*(pxc*uu(iy1,ix0)+px*uu(iy1,ix1))

  vq(i)=pyc*(pxc*vv(iy0,ix0)+px*vv(iy0,ix1)) &
      & +py*(pxc*vv(iy1,ix0)+px*vv(iy1,ix1))
enddo

return
end subroutine

!=======================================================================

subroutine con2grid(qq,xq,yq,dq,qavg,nextq,nptq,iopt)
! Contour -> grid conversion.  The contours are represented by
! nodes (xq(i),yq(i)), i = 1, ..., nptq, where nextq(i) gives the
! index of the node following i, dq is the jump in q across all
! contours, and qavg is average value of the field.

! The gridded field is returned in the array qq.

! The integer iopt is used for computing available potential energy (APE):
! iopt < 0: don't compute APE
! iopt = 0: initialise reference state potential energy and compute APE
! iopt > 0: compute APE (reference state PE assumed available in peref)

! The quantities peref and ape are in global commons

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(0:ny,0:nx)
double precision:: xq(npm),yq(npm)
integer:: nextq(npm)

 !Local arrays:
double precision:: qa(0:nyf+1,0:nxf)
double precision:: qjx(0:nxf+1),qbot(0:nxf)
double precision:: dx(nptq),dy(nptq)
double precision:: ybar(0:nlevm),area(nlevm)
integer:: ixc(nptq),nxc(nptq)
logical:: crossx(nptq)

if (nptq .eq. 0) then 
   !No contours to convert: return qq = 0:
  do ix=0,nx
    do iy=0,ny
      qq(iy,ix)=zero
    enddo
  enddo
  return
endif  

!------------------------------------------------------------------
 !Initialise interior x grid line crossing information and fill the
 !q jump array along lower boundary:
do i=1,nptq
  ixc(i)=int(one+dxxfi*(xq(i)-xbeg))
enddo
 !Here xbeg is very slightly larger than xmin so that a point on
 !the left edge has ixc = 0, but one with xq just greater than xbeg
 !has ixc = 1.  Similarly, a point on the right edge has ixc = nxf+1.
 !Note: dxxfi = dble(nxf)/((xmax-xmin)*(1-small)) where small = 1.d-12
 !(see casl.f90 initialisation)

do ix=0,nxf+1
  qjx(ix)=zero
enddo

do i=1,nptq
  ia=nextq(i)
  if (ia .gt. 0) then
     !A node with ia = 0 terminates a contour at a boundary
    dx(i)=xq(ia)-xq(i)
    dy(i)=yq(ia)-yq(i)
    nxc(i)=ixc(ia)-ixc(i)
    crossx(i)=(nxc(i) .ne. 0)
    if ((yq(ia)-ybeg)*(ybeg-yq(i)) .gt. zero) then
       !The contour segment (i,ia) crosses y = ybeg; find x location:
      py0=(ybeg-yq(i))/dy(i)
      ix=int(one+dxxfi*(xq(i)+py0*dx(i)-xbeg))
      qjx(ix)=qjx(ix)-dq*sign(one,dy(i))
       !Note: qjx gives the jump going from ix-1 to ix
    endif
  else
     !Here, there is no segment (i,next(i)) to consider:
    crossx(i)=.false.
  endif
enddo
 !Above, ybeg is very slightly greater than ymin to detect boundary crossings

 !Sum q jumps to obtain the gridded q along lower boundary:
qbot(0)=zero
 !Corner value cannot be determined a priori; qavg is used for this below
do ix=1,nxf
  qbot(ix)=qbot(ix-1)+qjx(ix)
enddo

!----------------------------------------------------------------
 !Initialise interior q jump array:
do ix=0,nxf
  do iy=0,nyf+1
    qa(iy,ix)=zero
  enddo
enddo

 !Determine x grid line crossings and accumulate q jumps:
do i=1,nptq
  if (crossx(i)) then
    jump=sign(1,nxc(i))
    ixbeg=ixc(i)+(jump-1)/2
    sdq=dq*sign(one,dx(i))
    ncr=0
    do while (ncr .ne. nxc(i)) 
      ix=ixbeg+ncr
      px0=(xxf(ix)-xq(i))/dx(i)
       !The contour crossed the fine grid line ix at the point
       !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
      iy=int(one+dyyfi*(yq(i)+px0*dy(i)-ybeg))
       !Increment q jump between the grid lines iy-1 & iy:
      qa(iy,ix)=qa(iy,ix)+sdq
       !Go on to consider next x grid line (if there is one):
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get q values by sweeping through y:
do ix=0,nxf
  qa(0,ix)=qbot(ix)
  do iy=1,nyf
    qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
  enddo
enddo

!------------------------------------------------------------------------
 !Possibly compute APE if iopt >= 0:
if (iopt .ge. 0) call getpe(qa,dq,qavg,iopt)

!------------------------------------------------------------------------

 !Average to inversion grid by repeated 1-2-1 averages in each direction:
call coarsen(qa,qq)
 !Restore average (qavg):
call average(qq,qavg0)

qadd=qavg-qavg0
do ix=0,nx
  do iy=0,ny
    qq(iy,ix)=qq(iy,ix)+qadd
  enddo
enddo
 !Now qq has the correct average

return
end subroutine

!=======================================================================

subroutine surgery(xq,yq,nextq,indq,npq,i1q,i2q,nq,nptq)
! Contour surgery.  The contours are represented by nodes (xq(i),yq(i)), 
! i = 1, ..., nptq, where nextq(i) gives the index of the node following 
! i, indq is an integer index distinguishing different contour levels, 
! npq is the number of points on a contour (returned only), i1q is the 
! starting index of a contour (returned only), i2q is the ending index
! (returned only), nq is the number of contours (returned only), and
! nptq is the total number of nodes (passed then modified and returned).

! After surgery, this routine redistributes nodes (renode).

! Major revision 01/01/2001 by D. G. Dritschel to accelerate surgery
! using local boxes to minimise search costs.  All surgery is now
! done in a single way.

! Now allows surgery with "open" contours which start and end in
! a boundary (completed 11 July 2012 by dgd @ St Andrews).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: xq(npm),yq(npm)
integer:: nextq(npm),indq(nm),npq(nm),i1q(nm),i2q(nm)

 !Local parameters and arrays:
integer,parameter:: ngbs=nx*ny
double precision:: xa(npm),ya(npm)
double precision:: xd(nprm),yd(nprm)
double precision:: dx(nptq),dy(nptq),dsq(nptq)
integer:: i1a(nm),i2a(nm)
integer:: jq1(nlevm),jq2(nlevm),iq1(nlevm),iq2(nlevm),levq(nlevm)
integer:: nspb(ngbs),kb1(ngbs),kb2(ngbs)
integer:: loc(nplm),list(nplm),node(nplm)
integer:: ibse(nq),ibss(nq),jbss(nq)
integer:: nexta(npm),icre(nm)

logical:: avail(nptq),edge(nptq)
logical:: open(nq),ave(nq),avs(nq)

if (nq .eq. 0) return

!--------------------------------------------------------------------
 !Calculate beginning and ending contours (jq1,jq2) for each distinct 
 !value of ind:
levp=indq(1)
nlev=1
jq1(1)=1
levq(1)=levp
do j=2,nq
  lev=indq(j)
  if (lev .gt. levp) then
    jq2(nlev)=j-1
    nlev=nlev+1
    jq1(nlev)=j
    levq(nlev)=lev
    levp=lev
  endif
enddo
jq2(nlev)=nq
 !Note: levq(lev) gives the q level (indq) of contours.

do lev=1,nlev
  iq1(lev)=i1q(jq1(lev))
  iq2(lev)=i2q(jq2(lev))
enddo

 !Get work arrays for efficient surgery:
do i=1,nptq-1
  dx(i)=xq(i+1)-xq(i)
  dy(i)=yq(i+1)-yq(i)
  dsq(i)=dx(i)**2+dy(i)**2
  avail(i)=.true.
  edge(i)=.false.
enddo
avail(nptq)=.true.
edge(nptq)=.false.

do j=1,nq
  ie=i2q(j)
  open(j)=(nextq(ie) .eq. 0)
  if (open(j)) then
     !Contour j is "open", i.e. it starts and ends at an edge
    avail(ie-2)=.false.
    avail(ie-1)=.false.
    avail(ie)  =.false.
    edge(ie-1)=.true.
    edge(ie)  =.true.
  else
     !Contour j is "closed"
    is=i1q(j)
    dx(ie)=xq(is)-xq(ie)
    dy(ie)=yq(is)-yq(ie)
    dsq(ie)=dx(ie)**2+dy(ie)**2
  endif
enddo
 !Note: avail(ib) is true when (i,nextq(i)), where i = nextq(ib), is 
 !a "proper" segment, i.e. one located fully within the domain, with
 !neither end attached to an edge.  Meanwhile, edge(ib) is false when
 !node i is not located at an edge.  When ib = ie, there is no node 
 !following ib when contour j is open.  The logicals avail and edge
 !are only used to separate and perform interior surgery first.
 !Thereafter, the meaning of the logical edge is changed - see below.

!-----------------------------------------------------------------
!-----------------------------------------------------------------

 !Begin a major loop over q levels:
nq=0
nptq=0

do lev=1,nlev
  ibeg=iq1(lev)
  iend=iq2(lev)

   !Work out optimal box number, nbox, for fast surgery:
  nptlev=iend-ibeg+1
   !nptlev: number of nodes having this q level.
  fnq=dble(nptlev)

   !Balance boxing costs, (3+15*nseg/nbox)*nbox, with surgery 
   !search costs, 24*nptlev*nseg/nbox, to work out optimal box number;
   !here we estimate nseg, the total number of segments counted
   !in all boxes, as (1+3*eps*nb)*nptlev, where nb=sqrt(nbox);
   !this has been verified in realistically complex tests).
   !This balance results in a quartic equation for a = nb/sqrt(nptlev)
   !whose solution only depends on r = 3*eps*sqrt(nptlev), where
   !eps=(average distance between adjacent nodes)/(domain area).
   !Fortunately, a only varies from 1.129 to 1.265 over the entire
   !range of r (from 0 to infinity).  Here, therefore, we simply
   !take a = 1.2, i.e. nb = 1.2*sqrt(nptlev), as an approximation:
  fnbx=1.2d0*sqrt(fnq*aspect)
   !aspect = x:y domain aspect ratio.
  fnby=fnbx/aspect
  nbx=max(min(nint(fnbx),nx),1)
   !nbx: number of boxes in the x direction; nbx <= nx.
  nby=max(min(nint(fnby),ny),1)
   !nby: number of boxes in the y direction; nby <= ny.
  nbox=nbx*nby
   !nbox: total number of boxes
  bwxi=oms*dble(nbx)/ellx
  bwyi=oms*dble(nby)/elly

   !Box all segments [i,nextq(i)] to reduce search costs in surgery:
  do mb=1,nbox
    nspb(mb)=0
  enddo
  nseg=0

  do ib=ibeg,iend
    if (avail(ib)) then
      i=nextq(ib)
      ia=nextq(i)
       !(i,ia) is a "proper" segment; find range of boxes spanned by it:
      if (dx(i) .gt. zero) then
        mbx1=int(bwxi*(xq(i) -xmin))
        mbx2=int(bwxi*(xq(ia)-xmin))
      else
        mbx1=int(bwxi*(xq(ia)-xmin))
        mbx2=int(bwxi*(xq(i) -xmin))
      endif
      if (dy(i) .gt. zero) then
        mby1=int(bwyi*(yq(i) -ymin))
        mby2=int(bwyi*(yq(ia)-ymin))
      else
        mby1=int(bwyi*(yq(ia)-ymin))
        mby2=int(bwyi*(yq(i) -ymin))
      endif
       !mbx1+1,mbx2+1 is the x range of boxes spanned by the segment while
       !mby1+1,mby2+1 is the y range.  

      do mbx=mbx1,mbx2
        do mby=mby1,mby2
          mb=nby*mbx+mby+1
          nspb(mb)=nspb(mb)+1
           !nspb counts number of segments that cross through box mb.
          nseg=nseg+1
           !nseg counts the total number of segments crossing all boxes.
          loc(nseg)=mb
           !loc gives the box location of the segment
          list(nseg)=ib
           !list gives the node before the one at the segment origin (i).
        enddo
      enddo
    endif
  enddo
   !nseg: the total number of segment crossings found for all boxes.
   !**NB: nseg/nbox = average number of segments crossing a box.

  kb1(1)=1
  do mb=1,nbox-1
    kb1(mb+1)=kb1(mb)+nspb(mb)
  enddo
  do mb=1,nbox
    kb2(mb)=kb1(mb)-1
  enddo

  do k=1,nseg
    mb=loc(k)
    ks=kb2(mb)+1
    node(ks)=list(k)
    kb2(mb)=ks
  enddo

   !With the above, segments [node(kb1(mb)),nextq(node(kb1(mb))],
   ![node(kb1(mb)+1),nextq(node(kb1(mb)+1)], ..., 
   ![node(kb2(mb)),nextq(node(kb2(mb))] cross through box mb.

   !Note: roughly (3+15*nseg/nbox)*nbox operations are required
   !      up to this point (in the loop over q levels) to do
   !      all the segment boxing.

!-------------------------------------------------------------------
   !Now search for surgery with segments crossing through the box 
   !containing i:
  do ib=ibeg,iend
    if (.not. edge(ib)) then
     !Node i = nextq(ib) is fully within the domain (not on an edge)
    i=nextq(ib)
    mb=nby*int(bwxi*(xq(i)-xmin))+int(bwyi*(yq(i)-ymin))+1
     !mb: box containing the node i.
    do k=kb1(mb),kb2(mb)
      isb=node(k)
      is=nextq(isb)
       !Segment [is,nextq(is)] lies in box mb.  Exclude segments 
       !ending in a edge or having node i at either endpoint:
      if (edge(is) .or. (is-ib)*(is-i) .eq. 0) goto 100
       !nextq(is) could lie at an edge as a result of previous surgery.

       !Next see if the node i lies within a distance dm to the line 
       !segment between is and nextq(is):
      delx=xq(is)-xq(i)
      dely=yq(is)-yq(i)
      aa=delx*dx(is)+dely*dy(is)

       !Note: roughly 24*nptlev*(nseg/nbox) operations are required 
       !up to the following statement (which is rarely satisfied); 
       !this is assumed to be the dominant cost of surgery.  
       !Even if each node i surgically reconnects on average 
       !once, the above estimate holds.

      if (aa*(aa+dsq(is))+d4small .lt. zero) then
         !Passing this condition, a perpendicular can be dropped 
         !onto the line segment from is to nextq(is).  
         ![Tests show that this is satisfied only 6.4% of the time.]
         !Check distance to line segment:
        cc=(delx*dy(is)-dely*dx(is))**2
         !NB: sqrt(cc/dsq(is)) is the distance to the segment and 
         !    dm2=dm**2 below.  This distance must be strictly 
         !    positive, yet less than dm, to permit surgery:
        if (cc*(cc-dm2*dsq(is)) .lt. zero) then
           !Surgery is now possible between node i and the segment 
           ![is,nextq(is)].
           ![Tests show that this is satisfied only 0.45% of the time.]
           !Move node i to a point midway between i and either is or 
           !nextq(is), whichever is closer:
          if (dsq(is)+two*aa .lt. zero) then
             !node i is closest to node nextq(is):
            dxa=f12*(delx+dx(is))
            dya=f12*(dely+dy(is))
            isb=is
            is=nextq(is)
          else
             !node i is closest to node is:
            dxa=f12*delx
            dya=f12*dely
          endif
           !Move node i & is to a common node; first deal with node i:
          xq(i)=xq(i)+dxa
          yq(i)=yq(i)+dya
          dx(i)=dx(i)-dxa
          dy(i)=dy(i)-dya
          dsq(i)=dx(i)**2+dy(i)**2
           !isb becomes the node before i:
          nextq(isb)=i
          dx(isb)=dx(isb)-dxa
          dy(isb)=dy(isb)-dya
          dsq(isb)=dx(isb)**2+dy(isb)**2

           !Now deal with node is:
          xq(is)=xq(i)
          yq(is)=yq(i)
          dx(is)=dx(is)+dxa
          dy(is)=dy(is)+dya
          dsq(is)=dx(is)**2+dy(is)**2
           !ib becomes the node before is:
          nextq(ib)=is
          dx(ib)=dx(ib)+dxa
          dy(ib)=dy(ib)+dya
          dsq(ib)=dx(ib)**2+dy(ib)**2

           !Now update i and look for further possible surgery:
          i=is
          if (i .eq. ib) goto 150
           !i = ib when a contour consists of a single node
        endif
      endif
      100 continue
    enddo
    150 continue
    endif
  enddo

!-------------------------------------------------------------------
   !Push nodes within a distance dm to a boundary to that boundary
   !and initialise logical array indicating that a node is available
   !for re-building contours below:
  do i=ibeg,iend
    if (abs(xmax-xq(i)) .lt. dm) then
      xq(i)=xmax
      edge(i)=.true.
    else if (abs(xq(i)-xmin) .lt. dm) then
      xq(i)=xmin
      edge(i)=.true.
    else if (abs(ymax-yq(i)) .lt. dm) then
      yq(i)=ymax
      edge(i)=.true.
    else if (abs(yq(i)-ymin) .lt. dm) then
      yq(i)=ymin
      edge(i)=.true.
    else
      edge(i)=.false.
    endif
    avail(i)=.true.
  enddo

   !Form lists of nodes (ibss,ibse) which start and end open contours:
  nse=0
  nss=0
  do ib=ibeg,iend
    i=nextq(ib)
    if (.not. edge(ib)) then
       !ib cannot end an open contour, and i > 0
      if (edge(i)) then
         !this node, i, lies on an edge of the domain
        ia=nextq(i)
        if (ia .gt. 0) then
           !there is a node after i
          if (edge(ia)) then
             !this node, ia, also lies on an edge; we thus stop the 
             !contour at i.  (ib,i) is a segment ending an open contour; 
             !record in a list:
            nse=nse+1
            ibse(nse)=ib
            ave(nse)=.true.
          endif
        else
           !there is no node after i.  (ib,i) is a segment ending an 
           !open contour; record in a list:
          nse=nse+1
          ibse(nse)=ib
          ave(nse)=.true.
        endif
      endif
    else
       !ib lies on the edge of the domain; it might end a contour, but
       !this possibility would be found above when the point before ib
       !is considered.
      if (i .gt. 0) then
         !there is a node after ib
        if (edge(i)) then
           !this node, i, also lies on an edge of the domain
          ia=nextq(i)
           !ensure ib is never used to start a contour:
          nextq(ib)=0
          avail(ib)=.false.
          if (ia .gt. 0) then
             !there is a node after i
            if (.not. edge(ia)) then
               !this node, ia, is in the interior; therefore (i,ia) is
               !a segment starting an open contour; record in a list:
              nss=nss+1
              ibss(nss)=i
              avs(nss)=.true.
            endif
          endif
        endif
      endif
    endif
  enddo

   !Finally check original open contours for segments starting a contour
   !not found above:
  do j=jq1(lev),jq2(lev)
    if (open(j)) then
       !contour j is open, and i below is its starting node (at an edge):
      i=i1q(j)
      if (avail(i)) then
         !This node is still available to start a contour; there must be
         !a node after (i.e. nextq(i)); therefore (i,nextq(i)) is a 
         !segment starting an open contour; record in a list:
        nss=nss+1
        ibss(nss)=i
        avs(nss)=.true.
      endif
       !ensure its original endpoint (unchanged by surgery above) is 
       !not available for starting a contour:
      avail(i2q(j))=.false.
    endif
  enddo
   !Note: there must be the same number of starting and ending segments,
   !      nss = nse.

   !Perform surgery between segments attached to a boundary:
  npe=0
  if (nse .gt. 0) then
     !There are open contours to process:
    do ke=1,nse
      if (ave(ke)) then
      ie=ibse(ke)
      iea=nextq(ie)
       !Do not allow the endpoint (iea) to ever start a contour:
      avail(iea)=.false.
      if (nextq(iea) .gt. 0) then
         !Make sure that any node following iea is not used to start
         !a *closed* contour (it could still start an open one):
        avail(nextq(iea))=.false.
        nextq(iea)=0
      endif

      do ks=1,nss
        if (avs(ks)) then
        is=ibss(ks)
        isa=nextq(is)

         !See if a perpendicular can be dropped from node is to the
         !segment (ie,iea):
        delx=xq(ie)-xq(is)
        dely=yq(ie)-yq(is)
        aa=delx*dx(ie)+dely*dy(ie)
        if (aa*(aa+dsq(ie))+d4small .lt. zero) then
           !Passing this condition, a perpendicular can be dropped 
           !onto the line segment (ie,iea).
          cc=(delx*dy(ie)-dely*dx(ie))**2
           !NB: sqrt(cc/dsq(ie)) is the distance to the segment and 
           !    dm2=dm**2 below.  This distance must be strictly 
           !    positive, yet less than dm, to permit surgery:
          if (cc*(cc-dm2*dsq(ie)) .lt. zero) then
             !Surgery is now possible between node is and the segment 
             !(ie,iea).  Move node is half-way along perpendicular and
             !eliminate node iea (avail(iea) = .false. above already):
            pe=-aa/dsq(ie)
            xq(is)=f12*(xq(is)+xq(ie)+pe*dx(ie))
            yq(is)=f12*(yq(is)+yq(ie)+pe*dy(ie))
            nextq(ie)=is

             !Don't consider these segments again in surgery:
            ave(ke)=.false.
            avs(ks)=.false.
            goto 200
          endif
        endif

         !See if a perpendicular can be dropped from node iea to the
         !segment (is,isa):
        delx=xq(is)-xq(iea)
        dely=yq(is)-yq(iea)
        aa=delx*dx(is)+dely*dy(is)
        if (aa*(aa+dsq(is))+d4small .lt. zero) then
           !Passing this condition, a perpendicular can be dropped 
           !onto the line segment (is,isa).
          cc=(delx*dy(is)-dely*dx(is))**2
           !NB: sqrt(cc/dsq(is)) is the distance to the segment and 
           !    dm2=dm**2 below.  This distance must be strictly 
           !    positive, yet less than dm, to permit surgery:
          if (cc*(cc-dm2*dsq(is)) .lt. zero) then
             !Surgery is now possible between node iea and the segment 
             !(is,isa).  Move node is half-way along perpendicular and
             !eliminate node iea (avail(iea) = .false. above already):
            ps=-aa/dsq(is)
            xq(is)=f12*(xq(iea)+xq(is)+ps*dx(is))
            yq(is)=f12*(yq(iea)+yq(is)+ps*dy(is))
            nextq(ie)=is

             !Don't consider these segments again in surgery:
            ave(ke)=.false.
            avs(ks)=.false.
            goto 200
          endif
        endif

        endif
      enddo
       !ends loop over ks

      endif
      200 continue
    enddo
     !ends loop over ke.

     !Next work out how many open contours remain by testing avs:
    do ks=1,nss
      if (avs(ks)) then
        npe=npe+1
        icre(npe)=ibss(ks)
      endif
    enddo
  endif

!-----------------------------------------------------------------------
   !It remains to rebuild the contours using the nextq() information

   !Current contour level:
  levt=levq(lev)

   !First deal with any open contours attached to boundaries:
  if (npe .gt. 0) then
    do ie=1,npe
       !A new contour (indexed nq) starts here:
      nq=nq+1

       !The starting node on the contour (coming out of a boundary):
      i=icre(ie)
      xd(1)=xq(i)
      yd(1)=yq(i)
      npd=1

      avail(i)=.false.
      is=nextq(i)
       !Generate the contour using nextq() and finish when the
       !other boundary is reached (is = 0):
      do while (is .ne. 0) 
        npd=npd+1
        xd(npd)=xq(is)
        yd(npd)=yq(is)
        avail(is)=.false.
        is=nextq(is)
      enddo

       !Re-distribute nodes on this contour:
      npq(nq)=npd
      if (npd .gt. 3) &
         & call renode_open(xd,yd,npd,xa(nptq+1),ya(nptq+1),npq(nq))
      if (npq(nq) .gt. 3) then
        i1a(nq)=nptq+1
        nptq=nptq+npq(nq)
        i2a(nq)=nptq
        indq(nq)=levt
        do is=i1a(nq),nptq-1
          nexta(is)=is+1
        enddo
        nexta(nptq)=0
      else
         !Delete contour if deemed too small (see renode):
        nq=nq-1
      endif

    enddo
  endif

   !Next deal with remaining closed contours:
  do i=ibeg,iend
    if (avail(i)) then
     !avail(i) is true if node i has not yet been associated with a contour.
     !start a new contour:
      nq=nq+1

     !First point on the contour:
      npd=1
      xd(1)=xq(i)
      yd(1)=yq(i)

      avail(i)=.false.
      is=nextq(i)
       !Generate the contour using nextq() and finish when the
       !original node (i) is reached:
      do while (is .ne. i) 
        npd=npd+1
        xd(npd)=xq(is)
        yd(npd)=yq(is)
        avail(is)=.false.
        is=nextq(is)
      enddo

       !Re-distribute nodes on this contour:
      npq(nq)=npd
      if (npd .gt. 3) &
         & call renode_closed(xd,yd,npd,xa(nptq+1),ya(nptq+1),npq(nq))
      if (npq(nq) .gt. 3) then
        i1a(nq)=nptq+1
        nptq=nptq+npq(nq)
        i2a(nq)=nptq
        indq(nq)=levt
        do is=i1a(nq),nptq-1
          nexta(is)=is+1
        enddo
        nexta(nptq)=i1a(nq)
      else
         !Delete contour if deemed too small (see renode):
        nq=nq-1
      endif

    endif
  enddo

   !Now go back and consider the next q level:
enddo
 !This ends the loop over q levels; surgery is complete.

 !Copy nodes back to those in the argument of the subroutine:
do i=1,nptq
  xq(i)=xa(i)
  yq(i)=ya(i)
  nextq(i)=nexta(i)
enddo
 !nextq(i) is not altered until surgery.

do j=1,nq
  i1q(j)=i1a(j)
  i2q(j)=i2a(j)
enddo

return
end subroutine

!==========================================================================

subroutine renode_closed(xd,yd,npd,xr,yr,npr)
! Re-nodes a single closed contour (xd(i),yd(i)), i = 1,...,npd 
!   and returns the new contour in (xr(i),yr(i)), i = 1,...,npr

! Note: a closed contour closes on itself

! If npr = 0, the contour is too small and should be removed

! -----------------------------------------------------------
! Intended for a rectangular domain with free-slip boundaries
! -----------------------------------------------------------

! ==> Constants zero, one, f12, f13, f16, small, small3, xmin, xmax,
! ==> ymin, ymax, elf, dmsq, dmi & densf are assumed defined

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: xd(nprm),yd(nprm),xr(nprm),yr(nprm)
 !Local parameters and arrays:
double precision:: dx(npd),dy(npd),a(npd),b(npd),c(npd),d(npd),e(npd)
double precision:: u(npd),v(npd)
integer:: node(npd+1),next(0:npd)
logical:: corner(npd)

!------------------------------------------------------------------
 !Define the node following a node:
do i=0,npd-1
  next(i)=i+1
enddo
next(npd)=1

 !Compute the contour increments:
do i=1,npd
  ia=next(i)
  dx(i)=xd(ia)-xd(i)
  dy(i)=yd(ia)-yd(i)
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npd
  v(i)=dx(i)*dx(i)+dy(i)*dy(i)
  e(i)=sqrt(v(i))
enddo

do ib=1,npd
  i=next(ib)
  u(i)=v(ib)
  a(i)=-dx(ib)
  c(i)=-dy(ib)
enddo

ncorn=0
do i=1,npd
  corner(i)=dx(i)*a(i)+dy(i)*c(i) .gt. zero
  if (corner(i)) then 
     !Keep track of corner locations for use in renoding below:
    ncorn=ncorn+1
    node(ncorn)=i
     !Set curvature to zero at corners:
    b(i)=zero
  else
    b(i)=(dx(i)*c(i)-a(i)*dy(i))/ &
     & sqrt((a(i)*v(i)-dx(i)*u(i))**2+(c(i)*v(i)-dy(i)*u(i))**2+small3)
  endif
enddo

do i=1,npd
  ia=next(i)
  u(i)=e(i)*(b(ia)+b(i))
  c(i)=e(i)*(b(ia)-b(i))
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npd
  a(i)=f16*c(i)-f12*u(i)
  b(i)=f12*(u(i)-c(i))
  c(i)=f13*c(i)
enddo

!------------------------------------------------------------------------
 !Use the spherical curvature expression (radius of the sphere = ell)
 !to ensure an adequate node density in low curvature regions.
do i=1,npd
  ww=one/(v(i)+dmsq)
  u(i)=ww*sqrt(elf*v(i)+u(i)**2)
  v(i)=ww*e(i)
enddo
 !NB: elf = 1/ell**2; v(i) = |xx_{i+1}-xx_{i}|**2; e(i)=sqrt{v(i)};
 !    u(i)/e(i) = (kappa_{i}+kappa_{i+1})/2; dmsq = (2*dm)**2

 !Re-assign curvature at a node from weighted average on either side
 !(v above is the weight):
do ib=1,npd
  i=next(ib)
  d(i)=(u(ib)+u(i))/(v(ib)+v(i)+small)
enddo

 !Re-average to get interval value (effectively, four curvature
 !values go into getting the final interval value, u(i)):
do i=1,npd
  ia=next(i)
  u(i)=f12*(d(i)+d(ia))
enddo

 !Compute fractional number of nodes to be placed between old
 !nodes i and i+1:
do i=1,npd
  e(i)=e(i)*min(dmi,densf*sqrt(u(i))+u(i))
enddo
 !NB: dmi = 2/delta; densf = 1/(amu*sqrt{ell})

sum=zero
do i=1,npd
  sum=sum+e(i)
enddo
 !Number of points on renoded contour:
npr=nint(sum)+1
if (npr .lt. 3) then
   !Contour too small - remove it:
  npr=0
  return
endif

!-------------------------------------------------------------
 !Redistribute nodes making sure to preserve corner locations:
if (ncorn .eq. 0) then
   !No corners - simplest case

   !Make the sum of e(i) equal to npr:
  fac=dble(npr)/sum
  do i=1,npd
    e(i)=fac*e(i)
  enddo

   !The first node is fixed; find the remaining node positions:
  xr(1)=xd(1)
  yr(1)=yd(1)
  acc=zero
  i=0
  do im=2,npr
    do while (acc .lt. one)
      i=i+1
      acc=acc+e(i)
    enddo
    acc=acc-one
    p=one-acc/e(i)
    eta=p*(a(i)+p*(b(i)+p*c(i)))
    xr(im)=xd(i)+p*dx(i)-eta*dy(i)
    yr(im)=yd(i)+p*dy(i)+eta*dx(i)
  enddo

else if (ncorn .eq. 1) then
   !A single corner - start new contour at the corner:

   !Make the sum of e(i) equal to npr:
  fac=dble(npr)/sum
  do i=1,npd
    e(i)=fac*e(i)
  enddo

   !The first node (the corner) is fixed; find the remaining nodes:
  i=node(1)
  xr(1)=xd(i)
  yr(1)=yd(i)
  acc=zero
  i=i-1
  do im=2,npr
    do while (acc .lt. one)
      i=next(i)
      acc=acc+e(i)
    enddo
    acc=acc-one
    p=one-acc/e(i)
    eta=p*(a(i)+p*(b(i)+p*c(i)))
    xr(im)=xd(i)+p*dx(i)-eta*dy(i)
    yr(im)=yd(i)+p*dy(i)+eta*dx(i)
  enddo

else
   !Multiple corners - start new contour at the first corner and
   !preserve locations of all corners:
  npr=0
  node(ncorn+1)=node(1)
  ibeg=node(1)

  do k=1,ncorn
    iend=node(k+1)
    i=ibeg
    sum=zero
    do while (i .ne. iend)
      sum=sum+e(i)
      i=next(i)
    enddo

    if (sum .gt. small) then
       !Adequate spacing exists between corners to renode this segment

       !Number of points on renoded contour segment:
      npseg=nint(sum)+1

       !Make the sum of e(i) equal to npseg
      fac=dble(npseg)/sum
      i=ibeg
      do while (i .ne. iend)
        e(i)=fac*e(i)
        i=next(i)
      enddo

       !Fix the first node along the segment as the corner:
      xr(npr+1)=xd(ibeg)
      yr(npr+1)=yd(ibeg)

      if (npseg .gt. 1) then
         !Find the remaining points along this segment:
        acc=zero
        i=ibeg-1
        do im=npr+2,npr+npseg
          do while (acc .lt. one)
            i=next(i)
            acc=acc+e(i)
          enddo
          acc=acc-one
          p=one-acc/e(i)
          eta=p*(a(i)+p*(b(i)+p*c(i)))
          xr(im)=xd(i)+p*dx(i)-eta*dy(i)
          yr(im)=yd(i)+p*dy(i)+eta*dx(i)
        enddo
      endif

      npr=npr+npseg

    endif
    ibeg=iend

     !go on and consider the next segment (if any) between corners:
  enddo
endif

 !Force points to remain inside the domain:
do i=1,npr
  xr(i)=min(xmax,max(xmin,xr(i)))
  yr(i)=min(ymax,max(ymin,yr(i)))
enddo

return
end subroutine

!=======================================================================

subroutine renode_open(xd,yd,npd,xr,yr,npr)
! Re-nodes a single open contour (xd(i),yd(i)), i = 1,...,npd 
! and returns the new contour in (xr(i),yr(i)), i = 1,...,npr

! Note: an open contour starts and ends in a boundary

! If npr = 0, the contour is too small and should be removed

! -----------------------------------------------------------
! Intended for a rectangular domain with free-slip boundaries
! -----------------------------------------------------------

! ==> Constants zero, one, f12, f13, f16, small, small3, xmin, xmax,
! ==> ymin, ymax, elf, dmsq, dmi & densf are assumed defined

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: xd(nprm),yd(nprm),xr(nprm),yr(nprm)
 !Local parameters and arrays:
double precision:: dx(npd),dy(npd),a(npd),b(npd),c(npd),d(npd),e(npd)
double precision:: u(npd),v(npd)
integer:: node(npd)
logical:: corner(npd)

!------------------------------------------------------------------
 !Compute the cubic interpolation coefficients:
do i=1,npd-1
  dx(i)=xd(i+1)-xd(i)
  a(i+1)=-dx(i)
  dy(i)=yd(i+1)-yd(i)
  c(i+1)=-dy(i)
  v(i)=dx(i)*dx(i)+dy(i)*dy(i)
  u(i+1)=v(i)
  e(i)=sqrt(v(i))
enddo

corner(1)=.true.
b(1)=zero
ncorn=1
node(1)=1
do i=2,npd-1
  corner(i)=dx(i)*a(i)+dy(i)*c(i) .gt. zero
  if (corner(i)) then 
     !Set curvature to zero at corners:
    b(i)=zero
     !Keep track of corner locations for use in renoding below:
    ncorn=ncorn+1
    node(ncorn)=i
  else
    b(i)=(dx(i)*c(i)-a(i)*dy(i))/ &
     & sqrt((a(i)*v(i)-dx(i)*u(i))**2+(c(i)*v(i)-dy(i)*u(i))**2+small3)
  endif
enddo
corner(npd)=.true.
b(npd)=zero
ncorn=ncorn+1
node(ncorn)=npd

do i=1,npd-1
  u(i)=e(i)*(b(i+1)+b(i))
  c(i)=e(i)*(b(i+1)-b(i))
  a(i)=f16*c(i)-f12*u(i)
  b(i)=f12*(u(i)-c(i))
  c(i)=f13*c(i)
enddo

!------------------------------------------------------------------------
 !Use the spherical curvature expression (radius of the sphere = ell)
 !to ensure an adequate node density in low curvature regions.
do i=1,npd-1
  ww=one/(v(i)+dmsq)
  u(i)=ww*sqrt(elf*v(i)+u(i)**2)
  v(i)=ww*e(i)
enddo
 !NB: elf = 1/ell**2; v(i) = |xx_{i+1}-xx_{i}|**2; e(i)=sqrt{v(i)};
 !    u(i)/e(i) = (kappa_{i}+kappa_{i+1})/2; dmsq = (2*dm)**2

 !Re-assign curvature at a node from weighted average on either side
 !(v above is the weight):
d(1)=u(1)/(v(1)+small)
do i=2,npd-1
  d(i)=(u(i-1)+u(i))/(v(i-1)+v(i)+small)
enddo
d(npd)=u(npd-1)/(v(npd-1)+small)
 !Note: endpoints only have one side to consider

 !Re-average to get interval value "curv" (effectively, four curvature
 !values go into getting this value), then compute the fractional number 
 !of nodes to be placed between old nodes i and i+1:
sum=zero
do i=1,npd-1
  curv=f12*(d(i)+d(i+1))
  e(i)=e(i)*min(dmi,densf*sqrt(curv)+curv)
  sum=sum+e(i)
enddo
 !NB: dmi = 2/delta; densf = 1/(amu*sqrt{ell})

 !Number of points on renoded contour:
npr=nint(sum)+1
if (npr .lt. 3) then
   !Contour too small - remove it:
  npr=0
  return
endif

!-------------------------------------------------------------
 !Redistribute nodes making sure to preserve corner locations:
npr=0
ibeg=1

do k=1,ncorn-1
  iend=node(k+1)
  sum=zero
  do i=ibeg,iend-1
    sum=sum+e(i)
  enddo

  if (sum .gt. small) then
     !Adequate spacing exists between corners to renode this segment

     !Number of points on renoded contour segment:
    npseg=nint(sum)+1

     !Make the sum of e(i) equal to npseg
    fac=dble(npseg)/sum
    do i=ibeg,iend-1
      e(i)=fac*e(i)
    enddo

     !Fix the first node along the segment as the corner:
    xr(npr+1)=xd(ibeg)
    yr(npr+1)=yd(ibeg)

    if (npseg .gt. 1) then
       !Find the remaining points along this segment:
      acc=zero
      i=ibeg-1
      do im=npr+2,npr+npseg
        do while (acc .lt. one)
          i=i+1
          acc=acc+e(i)
        enddo
        acc=acc-one
        p=one-acc/e(i)
        eta=p*(a(i)+p*(b(i)+p*c(i)))
        xr(im)=xd(i)+p*dx(i)-eta*dy(i)
        yr(im)=yd(i)+p*dy(i)+eta*dx(i)
      enddo
    endif

    npr=npr+npseg

  endif
  ibeg=iend

     !go on and consider the next segment (if any) between corners:
enddo

 !Last point:
npr=npr+1
xr(npr)=xd(npd)
yr(npr)=yd(npd)

 !Force points to remain inside the domain:
do i=1,npr
  xr(i)=min(xmax,max(xmin,xr(i)))
  yr(i)=min(ymax,max(ymin,yr(i)))
enddo

return
end subroutine

!======================================================================

subroutine getpe(qa,dq,qavg,iopt)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qa(0:nyf+1,0:nxf)

 !Local arrays:
double precision:: area(nlevm),sa(0:nlevm),ybar(0:nlevm)
double precision:: cygf(0:nyf),sxgf(0:nxf)
logical:: agtb

 !First get average qa on the fine grid:
qavg0=zero
do ix=1,nxf-1
  qavg0=qavg0+qa(0,ix)+qa(nyf,ix)
enddo
do iy=1,nyf-1
  qavg0=qavg0+qa(iy,0)+qa(iy,nxf)
enddo
qavg0=f12*qavg0+f14*(qa(0,0)+qa(nyf,0)+qa(0,nxf)+qa(nyf,nxf))
do ix=1,nxf-1
  do iy=1,nyf-1
    qavg0=qavg0+qa(iy,ix)
  enddo
enddo
qavg0=qavg0*dsumfi
 !dsumfi = 1/(nxf*nyf)
qadd=qavg-qavg0

if (iopt .eq. 0) then
   !Initialise reference state PE (peref):
  qamin=qa(0,0)
  qamax=qa(0,0)
  do ix=0,nxf
    do iy=0,nyf
      qamin=min(qamin,qa(iy,ix))
      qamax=max(qamax,qa(iy,ix))
    enddo
  enddo
  qoff=qadd+qamin

   !Compute the area occupied by each field level:
  dqi=one/dq
  nlev=nint((qamax-qamin)*dqi)
  do k=1,nlev+1
    area(k)=zero
  enddo

  k=nint((qa(0,0)-qamin)*dqi)+1
  area(k)=area(k)+f14
  k=nint((qa(0,nx)-qamin)*dqi)+1
  area(k)=area(k)+f14
  k=nint((qa(ny,0)-qamin)*dqi)+1
  area(k)=area(k)+f14
  k=nint((qa(ny,nx)-qamin)*dqi)+1
  area(k)=area(k)+f14

  do ix=1,nxf-1
    k=nint((qa(0,ix)-qamin)*dqi)+1
    area(k)=area(k)+f12
    k=nint((qa(ny,ix)-qamin)*dqi)+1
    area(k)=area(k)+f12
  enddo

  do iy=1,nyf-1
    k=nint((qa(iy,0)-qamin)*dqi)+1
    area(k)=area(k)+f12
    k=nint((qa(iy,nx)-qamin)*dqi)+1
    area(k)=area(k)+f12
  enddo

  do ix=1,nxf-1
    do iy=1,nyf-1
      k=nint((qa(iy,ix)-qamin)*dqi)+1
      area(k)=area(k)+one
    enddo
  enddo

  if (.not. tilt) then 
     !The domain is level.
     !Compute the height of each field level by summing areas / L_x:
    yfac=gareaf/ellx
    ybar(0)=zero
    do k=1,nlev
      ybar(k)=ybar(k-1)+yfac*area(k)
    enddo
    ybar(nlev+1)=elly

     !Adjust the field so it has the correct average and compute peref:
    peref=zero
    do k=0,nlev
      peref=peref+(qoff+dble(k)*dq)*(ybar(k+1)**2-ybar(k)**2)
    enddo
    peref=qavg*domarea*hly-peref*hlx
  else
     !Domain is tilted in this case:
    astheta=abs(stheta)
    s2theta=two*astheta*ctheta
    csc2theta=one/s2theta

    a=ellx*astheta
    b=elly*ctheta
    agtb=(a .gt. b)

    if (agtb) then
      acor=f12*elly**2*ctheta/astheta
      yclo=b
      ychi=a
      wmid=elly/astheta
    else
      acor=f12*ellx**2*astheta/ctheta
      yclo=a
      ychi=b
      wmid=ellx/ctheta
    endif
    amid=domarea-two*acor
    yfac=one/wmid

    sa(0)=zero
    do k=1,nlev
      sa(k)=sa(k-1)+gareaf*area(k)
    enddo
    sa(nlev+1)=domarea

    ybar(0)=zero
    do k=1,nlev
      if (sa(k) .lt. acor) then
        ybar(k)=sqrt(sa(k)*s2theta)
      else if (sa(k) .lt. acor+amid) then
        ybar(k)=yclo+(sa(k)-acor)*yfac
      else
        ybar(k)=a+b-sqrt((domarea-sa(k))*s2theta)
      endif
    enddo
    ybar(nlev+1)=a+b

     !Compute reference state potential energy by explicit integration 
     !over the tilted domain:
    peref=zero
    do k=1,nlev+1
      qlev=qoff+dq*dble(k-1)
      if (sa(k) .lt. acor) then
        peref=peref+qlev*f23*csc2theta*(ybar(k)**3-ybar(k-1)**3)
      else if (sa(k) .lt. acor+amid) then
        if (sa(k-1) .ge. acor) then
          peref=peref+qlev*f12*wmid*(ybar(k)**2-ybar(k-1)**2)
        else
          peref=peref+qlev*(f23*csc2theta*(yclo**3-ybar(k-1)**3)+ &
                          & f12*wmid*(ybar(k)**2-yclo**2))
        endif
      else
        if (sa(k-1) .ge. acor+amid) then
          peref=peref+qlev*csc2theta*((a+b)*(ybar(k)**2-ybar(k-1)**2)- &
                                    & f23*(ybar(k)**3-ybar(k-1)**3))
        else if (sa(k-1) .lt. acor) then
          peref=peref+qlev*(f23*csc2theta*(yclo**3-ybar(k-1)**3)+ &
                            f12*wmid*(ychi**2-yclo**2)+ &
                          & csc2theta*((a+b)*(ybar(k)**2-ychi**2)- &
                          & f23*(ybar(k)**3-ychi**3)))
        else
          peref=peref+qlev*(f12*wmid*(ychi**2-ybar(k-1)**2)+ &
                          & csc2theta*((a+b)*(ybar(k)**2-ychi**2)- &
                          & f23*(ybar(k)**3-ychi**3)))
        endif
      endif
    enddo
    peref=qavg*domarea*(hly*ctheta+hlx*astheta)-peref
  endif
endif

if (.not. tilt) then
   !Compute PE by integrating (qa+qadd)*(ymin-y) over the domain:
  pe=zero
  do ix=1,nxf-1
    pe=pe+(qa(nyf,ix)+qadd)*elly
  enddo
  do iy=1,nyf-1
    pe=pe+(qa(iy,0)+qadd+qa(iy,nxf)+qadd)*ygf(iy)
  enddo
  pe=f12*pe+f14*(qa(nyf,0)+qadd+qa(nyf,nxf)+qadd)*elly
  do ix=1,nxf-1
    do iy=1,nyf-1
      pe=pe+(qa(iy,ix)+qadd)*ygf(iy)
    enddo
  enddo
  pe=-pe*gareaf
else 
   !Compute PE by integrating (qa+qadd)*[(ymin-y)*cos(theta)+(xmin-x)*sin(theta)]
   ! over the domain:
  do ix=0,nxf
    sxgf(ix)=stheta*xgf(ix)
  enddo
  do iy=0,nyf
    cygf(iy)=ctheta*ygf(iy)
  enddo

  pe=zero
  do ix=1,nxf-1
    pe=pe+(qa(nyf,ix)+qadd)*(cygf(nyf)+sxgf(ix)) &
       & +(qa(0,ix)  +qadd)*(cygf(0)  +sxgf(ix)) 
  enddo
  do iy=1,nyf-1
    pe=pe+(qa(iy,0)  +qadd)*(cygf(iy)+sxgf(0)) &
       & +(qa(iy,nxf)+qadd)*(cygf(iy)+sxgf(nxf))
  enddo
  pe=f12*pe+f14*((qa(0,0)    +qadd)*(cygf(0)  +sxgf(0))   &
              & +(qa(nyf,0)  +qadd)*(cygf(nyf)+sxgf(0))   &
              & +(qa(0,nxf)  +qadd)*(cygf(0)  +sxgf(nxf)) &
              & +(qa(nyf,nxf)+qadd)*(cygf(nyf)+sxgf(nxf))) 
  do ix=1,nxf-1
    do iy=1,nyf-1
      pe=pe+(qa(iy,ix)+qadd)*(cygf(iy)+sxgf(ix))
    enddo
  enddo
  pe=-pe*gareaf
endif

if (osci) then 
   !Define potential energy relative to the domain centre:
  ape=pe+qavg*domarea*(hly*ctheta+hlx*stheta)
else
   !Get APE by subtracting reference state PE:
  ape=pe-peref
endif
  
return
end subroutine

!======================================================================

subroutine coarsen(qa,qq)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qa(0:nyf+1,0:nxf),qq(0:ny,0:nx)

 !Local array:
double precision:: qh(0:nyf/2,0:nxf/2)

 !***Warning: The coarsening assumes nxf = 4*nx, nyf = 4*ny
nxh=nxf/2
nyh=nyf/2

do ix=0,nxh
  mix=2*ix
  mixp1=nxf-abs(nxf-abs(mix+1))
  mixm1=nxf-abs(nxf-abs(mix-1))
  do iy=0,nyh
    miy=2*iy
    miyp1=nyf-abs(nyf-abs(miy+1))
    miym1=nyf-abs(nyf-abs(miy-1))
    qh(iy,ix)=f14*qa(miy,mix) &
         &  +f18*(qa(miy,mixm1)+qa(miy,mixp1)  &
         &       +qa(miym1,mix)+qa(miyp1,mix)) &
         & +f116*(qa(miym1,mixm1)+qa(miym1,mixp1)  &
         &       +qa(miyp1,mixm1)+qa(miyp1,mixp1))
  enddo
enddo

do ix=0,nx
  mix=2*ix
  mixp1=nxh-abs(nxh-abs(mix+1))
  mixm1=nxh-abs(nxh-abs(mix-1))
  do iy=0,ny
    miy=2*iy
    miyp1=nyh-abs(nyh-abs(miy+1))
    miym1=nyh-abs(nyh-abs(miy-1))
    qq(iy,ix)=f14*qh(miy,mix) &
         &  +f18*(qh(miy,mixm1)+qh(miy,mixp1)  &
         &       +qh(miym1,mix)+qh(miyp1,mix)) &
         & +f116*(qh(miym1,mixm1)+qh(miym1,mixp1)  &
         &       +qh(miyp1,mixm1)+qh(miyp1,mixp1))
  enddo
enddo

return
end subroutine

!======================================================================

 !Main end module
end module
