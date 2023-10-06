module contours

! Module contains all subroutines related to contour advection, for the 
! casl suite of f90 codes in a doubly-periodic domain

use constants

implicit none

 !PV contours:
double precision:: xq(npm),yq(npm),qjump
integer:: indq(nm),npq(nm),i1q(nm),i2q(nm)
integer:: nextq(npm),nptq,nq,npta,na

 !Contour to grid conversion quantities (fine grid used in this module):
double precision:: xgf(nxf),ygf(nyf)
double precision,parameter:: glxf=ellx/dble(nxf),glxfi=dble(nxf)/ellx
double precision,parameter:: glyf=elly/dble(nyf),glyfi=dble(nyf)/elly

 !Contour to grid conversion quantities (ultra-fine grid used in congen):
double precision:: xgu(nxu),ygu(nyu)
double precision,parameter:: glxu=ellx/dble(nxu),glxui=dble(nxu)/ellx
double precision,parameter:: glyu=elly/dble(nyu),glyui=dble(nyu)/elly

 !Next grid points used in bilinear interpolation:
integer:: ixp(nx),iyp(ny)

 !Area weights used for interpolation in module congen:
double precision:: w00(mgu,mgu),w10(mgu,mgu),w01(mgu,mgu),w11(mgu,mgu)
integer:: ixfw(nxu),ix0w(nxu),ix1w(nxu)
integer:: iyfw(nyu),iy0w(nyu),iy1w(nyu)

 !Grid box reference indices used in congen:
integer:: ibx(nxu+1),iby(nyu+1)

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

!------------------------------------------------------------------------
 !Next grid points used in velocity interpolation (velint) and elsewhere; 
 !enforces periodicity in x and y:
do ix=1,nx-1
  ixp(ix)=ix+1
enddo
ixp(nx)=1

do iy=1,ny-1
  iyp(iy)=iy+1
enddo
iyp(ny)=1

 !Fine x-grid lines needed for contour-to-grid conversion (con2grid):
do ix=1,nxf
  xgf(ix)=xmin+glxf*dble(ix-1)
enddo

do iy=1,nyf
  ygf(iy)=ymin+glyf*dble(iy-1)
enddo

!==========================================================================
 !Initialise ultra-fine grid lines needed in contour regeneration (congen):
do ix=1,nxu
  xgu(ix)=xmin+glxu*dble(ix-1)
enddo

do iy=1,nyu
  ygu(iy)=ymin+glyu*dble(iy-1)
enddo

 !Grid box reference indices used in congen:
ibx(1)=nxu-1
do ix=2,nxu+1
  ibx(ix)=ix-2
enddo
iby(1)=nyu-1
do iy=2,nyu+1
  iby(iy)=iy-2
enddo

!----------------------------------------------------------------------
 !Area weights for interpolation of a gridded field onto the ultra-fine 
 !horizontal grid (congen):
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

do ix=1,nxu
  ixx=(ix-1)/mgu
  ix0w(ix)=1+ixx
  ix1w(ix)=2+ixx-nx*(ix0w(ix)/nx)
  ixfw(ix)=ix-mgu*ixx
enddo

do iy=1,nyu
  iyy=(iy-1)/mgu
  iy0w(iy)=1+iyy
  iy1w(iy)=2+iyy-ny*(iy0w(iy)/ny)
  iyfw(iy)=iy-mgu*iyy
enddo

return
end subroutine

!=======================================================================
      
subroutine velint(uu,vv,uq,vq)

! Bi-linearly interpolates current velocity (uu,vv) to the contour
! nodes (xq,yq) and stores the result in (uq,vq).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: uu(ny,nx),vv(ny,nx)
double precision:: uq(nptq),vq(nptq)

do i=1,nptq
  xx=glxi*(xq(i)-xmin)
  ix0=1+int(xx)
  pxc=dble(ix0)-xx
  px=one-pxc
  ix1=ixp(ix0)

  yy=glyi*(yq(i)-ymin)
  iy0=1+int(yy)
  pyc=dble(iy0)-yy
  py=one-pyc
  iy1=iyp(iy0)

  uq(i)=pyc*(pxc*uu(iy0,ix0)+px*uu(iy0,ix1)) &
        +py*(pxc*uu(iy1,ix0)+px*uu(iy1,ix1))

  vq(i)=pyc*(pxc*vv(iy0,ix0)+px*vv(iy0,ix1)) &
        +py*(pxc*vv(iy1,ix0)+px*vv(iy1,ix1))
enddo

return
end subroutine

!=======================================================================

subroutine con2grid(qc)
! Contour -> grid conversion.
! The gridded field is returned in the array qc.
! Note: the PV anomaly (PV - beta*y) is returned

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qc(ny,nx)

 !Local arrays:
double precision:: qod0(nyf/2),qod1(nyf/2),qod2(nyf/2)
double precision:: qev0(0:nyf/2),qev1(0:nyf/2),qev2(0:nyf/2)
double precision:: qa(nyf+1,nxf)
double precision:: qjx(nxf),qbot(nxf)
double precision:: dx(npta),dy(npta)
integer:: ixc(npta),nxc(npta)
logical:: crossx(npta)

if (npta .eq. 0) then 
   !No contours to convert: return qc = 0:
  qc=zero
  return
endif  

!------------------------------------------------------------------
 !Initialise interior x grid line crossing information and fill the
 !q jump array along lower boundary:
do i=1,npta
  ixc(i)=1+int(glxfi*(xq(i)-xmin))
enddo

qjx=zero

do i=1,npta
  ia=nextq(i)
  xx=xq(ia)-xq(i)
  dx(i)=xx-ellx*dble(int(xx*hlxi))
  yy=yq(ia)-yq(i)
  dy(i)=yy-elly*dble(int(yy*hlyi))
  ixdif=ixc(ia)-ixc(i)
  nxc(i)=ixdif-nxf*((2*ixdif)/nxf)
  crossx(i)=(nxc(i) .ne. 0)
  if (abs(yy) .gt. hly) then
     !The contour segment (i,ia) crosses y = ymin; find x location:
    cc=sign(one,dy(i))
    p=-(yq(i)-hly*cc)/(dy(i)+small)
    ix=1+int(glxfi*(mod(ellx+hlx+xq(i)+p*dx(i),ellx)))
    qjx(ix)=qjx(ix)-qjump*cc
     !Note: qjx gives the jump going from ix to ix+1
  endif
enddo

 !Sum q jumps to obtain the gridded q along lower boundary:
qbot(1)=zero
 !Corner value cannot be determined a priori; qavg is used for this below
do ix=1,nxf-1
  qbot(ix+1)=qbot(ix)+qjx(ix)
enddo

!----------------------------------------------------------------
 !Initialise interior q jump array:
qa=zero

 !Determine x grid line crossings and accumulate q jumps:
do i=1,npta
  if (crossx(i)) then
    jump=sign(1,nxc(i))
    ixbeg=ixc(i)+(jump-1)/2+nxf
    sqjump=qjump*sign(one,dx(i))
    ncr=0
    do while (ncr .ne. nxc(i)) 
      ix=1+mod(ixbeg+ncr,nxf)
      xx=xgf(ix)-xq(i)
      px0=(xx-ellx*dble(int(xx*hlxi)))/dx(i)
       !The contour crossed the fine grid line ix at the point
       !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
      yy=yq(i)+px0*dy(i)
      iy=2+int(glyfi*(yy-elly*dble(int(yy*hlyi))-ymin))
       !Increment q jump between the grid lines iy-1 & iy:
      qa(iy,ix)=qa(iy,ix)+sqjump
       !Go on to consider next x grid line (if there is one):
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get q values by sweeping through y:
do ix=1,nxf 
  qa(1,ix)=qbot(ix)
  do iy=2,nyf
    qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
  enddo
enddo

if (beffect) then
   !Remove beta*y to define anomaly before averaging below:
  do ix=1,nxf 
    do iy=1,nyf
      qa(iy,ix)=qa(iy,ix)-beta*ygf(iy)
    enddo
  enddo
endif

!------------------------------------------------------------------------
 !Average the PV field in qa to the coarser grid (ny,nx):
nxh=nxf
nyh=nyf
do while (nxh .gt. nx)
  nxff=nxh
  nxh=nxh/2
  nyff=nyh
  nyh=nyh/2
   !Perform nine-point averaging:
  do iy=1,nyh
    miy=2*iy
    qod2(iy)=qa(miy-1,nyff)
    qev2(iy)=qa(miy,nyff)
  enddo
  qev2(0)=qa(nyff,nyff)
  do ix=1,nxh
    mix=2*ix
    mixm1=mix-1
    do iy=1,nyh
      miy=2*iy
      qod1(iy)=qa(miy-1,mixm1)
      qod0(iy)=qa(miy-1,mix)
      qev1(iy)=qa(miy,mixm1)
      qev0(iy)=qa(miy,mix)
    enddo
    qev1(0)=qev1(nyh)
    qev0(0)=qev0(nyh)
    do iy=1,nyh
      qa(iy,ix)=0.0625d0*(qev0(iy)+qev0(iy-1)+qev2(iy)+qev2(iy-1)) &
              & +0.125d0*(qev1(iy)+qev1(iy-1)+qod0(iy)+qod2(iy)) &
              &   +0.25d0*qod1(iy)
    enddo
    do iy=1,nyh
      qod2(iy)=qod0(iy)
      qev2(iy)=qev0(iy)
    enddo
    qev2(0)=qev0(0)
  enddo
enddo

 !Remove average PV and copy into qc array:
qavg0=dsumi*sum(qa(1:ny,1:nx))
qc=qa(1:ny,1:nx)-qavg0
 !Now qc has zero average

return
end subroutine

!=======================================================================

subroutine surgery
! Contour surgery and node redistribution.
! Revised from the channel geometry stratified case 30 Sep 2013. 

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local parameters and arrays:
integer,parameter:: ngbs=nx*ny
double precision:: xa(npm),ya(npm)
double precision:: xd(nprm),yd(nprm)
double precision:: dx(nptq),dy(nptq),dsq(nptq)
integer:: i1a(nm),i2a(nm)
integer:: jq1(nlevm),jq2(nlevm),iq1(nlevm),iq2(nlevm),levq(nlevm)
integer:: nspb(ngbs),kb1(ngbs),kb2(ngbs)
integer:: loc(nplm),list(nplm),node(nplm)
integer:: nexta(npm),icre(nm)
logical:: avail(nptq)

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
  xx=xq(i+1)-xq(i)
  dx(i)=xx-ellx*dble(int(xx*hlxi))
  yy=yq(i+1)-yq(i)
  dy(i)=yy-elly*dble(int(yy*hlyi))
  dsq(i)=dx(i)**2+dy(i)**2
  avail(i)=.true.
enddo
avail(nptq)=.true.

do j=1,nq
  ie=i2q(j)
  is=i1q(j)
  xx=xq(is)-xq(ie)
  dx(ie)=xx-ellx*dble(int(xx*hlxi))
  yy=yq(is)-yq(ie)
  dy(ie)=yy-elly*dble(int(yy*hlyi))
  dsq(ie)=dx(ie)**2+dy(ie)**2
enddo

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

  dnbx=dble(nbx)
  dnby=dble(nby)
  do ib=ibeg,iend
    i=nextq(ib)
     !Find range of boxes spanned by the segment (i,nextq(i)):
    if (dx(i) .gt. zero) then
      mbx1=int(dnbx+bwxi*(xq(i)-xmin))
      mbx2=int(dnbx+bwxi*(xq(i)+dx(i)-xmin))
    else
      mbx1=int(dnbx+bwxi*(xq(i)+dx(i)-xmin))
      mbx2=int(dnbx+bwxi*(xq(i)-xmin))
    endif
    if (dy(i) .gt. zero) then
      mby1=int(dnby+bwyi*(yq(i)-ymin))
      mby2=int(dnby+bwyi*(yq(i)+dy(i)-ymin))
    else
      mby1=int(dnby+bwyi*(yq(i)+dy(i)-ymin))
      mby2=int(dnby+bwyi*(yq(i)-ymin))
    endif
     !mbx1+1,mbx2+1 is the x range of boxes spanned by the segment while
     !mby1+1,mby2+1 is the y range.  

    do mbx=mbx1,mbx2
      mmbx=mod(mbx,nbx)
      do mby=mby1,mby2
        mb=nby*mmbx+mod(mby,nby)+1
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
    i=nextq(ib)
    mb=nby*int(bwxi*(xq(i)-xmin))+int(bwyi*(yq(i)-ymin))+1
     !mb: box containing the node i.
    do k=kb1(mb),kb2(mb)
      isb=node(k)
      is=nextq(isb)
       !Segment [is,nextq(is)] lies in box mb.  Exclude segments 
       !ending in a edge or having node i at either endpoint:
      if ((is-ib)*(is-i) .eq. 0) cycle
       !nextq(is) could lie at an edge as a result of previous surgery.

       !Next see if the node i lies within a distance dm to the line 
       !segment between is and nextq(is):
      xx=xq(is)-xq(i)
      delx=xx-ellx*dble(int(xx*hlxi))
      yy=yq(is)-yq(i)
      dely=yy-elly*dble(int(yy*hlyi))
      dotp=delx*dx(is)+dely*dy(is)

       !Note: roughly 24*nptlev*(nseg/nbox) operations are required 
       !up to the following statement (which is rarely satisfied); 
       !this is assumed to be the dominant cost of surgery.  
       !Even if each node i surgically reconnects on average 
       !once, the above estimate holds.

      if (dotp*(dotp+dsq(is))+d4small .lt. zero) then
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
          if (dsq(is)+two*dotp .lt. zero) then
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
          xx=xq(i)+dxa
          xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
          yy=yq(i)+dya
          yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
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
          if (i .eq. ib) exit
           !i = ib when a contour consists of a single node
        endif
      endif
    enddo
  enddo

!-----------------------------------------------------------------------
   !It remains to rebuild the contours using the nextq() information

   !Current contour level:
  levt=levq(lev)

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
         & call renode(xd,yd,npd,xa(nptq+1),ya(nptq+1),npq(nq))
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

 !Find number of active PV contours and nodes:
if (tracer) then
   !There are tracer contours with indq(j) = 9999:
  j=nq
  do while (indq(j) .eq. 9999)
    j=j-1
  enddo
  na=j
  npta=i2q(j)
else
   !No tracer contours: all are "active":
  na=nq
  npta=nptq
endif

return
end subroutine

!==========================================================================

subroutine renode(xd,yd,npd,xr,yr,npr)
! Re-nodes a single closed contour (xd(i),yd(i)), i = 1,...,npd 
!   and returns the new contour in (xr(i),yr(i)), i = 1,...,npr

! Note: a closed contour closes on itself

! If npr = 0, the contour is too small and should be removed

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
  xx=xd(ia)-xd(i)
  dx(i)=xx-ellx*dble(int(xx*hlxi))
  yy=yd(ia)-yd(i)
  dy(i)=yy-elly*dble(int(yy*hlyi))
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
    disn=p*(a(i)+p*(b(i)+p*c(i)))
    xx=xd(i)+p*dx(i)-disn*dy(i)
    xr(im)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yy=yd(i)+p*dy(i)+disn*dx(i)
    yr(im)=oms*(yy-elly*dble(int(yy*hlyi)))
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
    disn=p*(a(i)+p*(b(i)+p*c(i)))
    xx=xd(i)+p*dx(i)-disn*dy(i)
    xr(im)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yy=yd(i)+p*dy(i)+disn*dx(i)
    yr(im)=oms*(yy-elly*dble(int(yy*hlyi)))
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
          disn=p*(a(i)+p*(b(i)+p*c(i)))
          xx=xd(i)+p*dx(i)-disn*dy(i)
          xr(im)=oms*(xx-ellx*dble(int(xx*hlxi)))
          yy=yd(i)+p*dy(i)+disn*dx(i)
          yr(im)=oms*(yy-elly*dble(int(yy*hlyi)))
        enddo
      endif

      npr=npr+npseg

    endif
    ibeg=iend

     !go on and consider the next segment (if any) between corners:
  enddo
endif

return
end subroutine

!=======================================================================

subroutine circulation(uu,vv,ax,ay,jj,circh,circm,dchdt)

! Computes circulations C_h = int{omega dxdy} = int_C{u*dx}
! and C_m = int{j dxdy} = int_C{B*dx}, together with dC_h/dt

implicit none

 !Passed variables:
double precision:: uu(ny,nx),vv(ny,nx),ax(ny,nx),ay(ny,nx),jj(ny,nx)
double precision:: circh,circm,dchdt

 !Local variables:
double precision:: dx(npta+1:nptq),dy(npta+1:nptq)
double precision::  b(npta+1:nptq), c(npta+1:nptq)
double precision:: spl,dspl,xx,pxc,px,yy,pyc,py
double precision:: x1,y1,dx1,dy1,uu1,vv1,ax1,ay1,jj1
double precision:: x2,y2,dx2,dy2,uu2,vv2,ax2,ay2,jj2
integer:: ibeg,iend,i,ix0,ix1,iy0,iy1

!----------------------------------------------------------------
 !First calculate the cubic interpolation coefficients:
ibeg=npta+1
iend=nptq
call cubic(ibeg,iend,dx,dy,b,c)

 !Next accumulate contour integrals:
circh=zero
dchdt=zero
circm=zero
do i=ibeg,iend
   !1st Gaussian point:
  spl=gp1*(gp1*(b(i)+gp1*c(i))-b(i)-c(i))
  xx=xq(i)+gp1*dx(i)-spl*dy(i)
  x1=xx-ellx*dble(int(xx*hlxi))
  yy=yq(i)+gp1*dy(i)+spl*dx(i)
  y1=yy-elly*dble(int(yy*hlyi))
  dspl=gp1*(two*b(i)+three*gp1*c(i))-b(i)-c(i)
  dx1=dx(i)-dspl*dy(i)
  dy1=dy(i)+dspl*dx(i)
   !Interpolate u, v, A_x, A_y & Lap(A) = -j:
  xx=glxi*(x1-xmin)
  ix0=1+int(xx)
  pxc=dble(ix0)-xx
  px=one-pxc
  ix1=ixp(ix0)

  yy=glyi*(y1-ymin)
  iy0=1+int(yy)
  pyc=dble(iy0)-yy
  py=one-pyc
  iy1=iyp(iy0)

  uu1=pyc*(pxc*uu(iy0,ix0)+px*uu(iy0,ix1))+py*(pxc*uu(iy1,ix0)+px*uu(iy1,ix1))
  vv1=pyc*(pxc*vv(iy0,ix0)+px*vv(iy0,ix1))+py*(pxc*vv(iy1,ix0)+px*vv(iy1,ix1))
  ax1=pyc*(pxc*ax(iy0,ix0)+px*ax(iy0,ix1))+py*(pxc*ax(iy1,ix0)+px*ax(iy1,ix1))
  ay1=pyc*(pxc*ay(iy0,ix0)+px*ay(iy0,ix1))+py*(pxc*ay(iy1,ix0)+px*ay(iy1,ix1))
  jj1=pyc*(pxc*jj(iy0,ix0)+px*jj(iy0,ix1))+py*(pxc*jj(iy1,ix0)+px*jj(iy1,ix1))

   !2nd Gaussian point:
  spl=gp2*(gp2*(b(i)+gp2*c(i))-b(i)-c(i))
  xx=xq(i)+gp2*dx(i)-spl*dy(i)
  x2=xx-ellx*dble(int(xx*hlxi))
  yy=yq(i)+gp2*dy(i)+spl*dx(i)
  y2=yy-elly*dble(int(yy*hlyi))
  dspl=gp2*(two*b(i)+three*gp2*c(i))-b(i)-c(i)
  dx2=dx(i)-dspl*dy(i)
  dy2=dy(i)+dspl*dx(i)
   !Interpolate u, v, A_x, A_y & Lap(A) = -j:
  xx=glxi*(x2-xmin)
  ix0=1+int(xx)
  pxc=dble(ix0)-xx
  px=one-pxc
  ix1=ixp(ix0)

  yy=glyi*(y2-ymin)
  iy0=1+int(yy)
  pyc=dble(iy0)-yy
  py=one-pyc
  iy1=iyp(iy0)

  uu2=pyc*(pxc*uu(iy0,ix0)+px*uu(iy0,ix1))+py*(pxc*uu(iy1,ix0)+px*uu(iy1,ix1))
  vv2=pyc*(pxc*vv(iy0,ix0)+px*vv(iy0,ix1))+py*(pxc*vv(iy1,ix0)+px*vv(iy1,ix1))
  ax2=pyc*(pxc*ax(iy0,ix0)+px*ax(iy0,ix1))+py*(pxc*ax(iy1,ix0)+px*ax(iy1,ix1))
  ay2=pyc*(pxc*ay(iy0,ix0)+px*ay(iy0,ix1))+py*(pxc*ay(iy1,ix0)+px*ay(iy1,ix1))
  jj2=pyc*(pxc*jj(iy0,ix0)+px*jj(iy0,ix1))+py*(pxc*jj(iy1,ix0)+px*jj(iy1,Ix1))

   !Accumulate hydrodynamic circulation:
  circh=circh+uu1*dx1+vv1*dy1+uu2*dx2+vv2*dy2
   !Accumulate magnetic circulation:
  circm=circm-ax1*dy1+ay1*dx1-ax2*dy2+ay2*dx2
   !Accumulate hydrodynamic circulation rate of change:
  dchdt=dchdt+jj1*(ax1*dx1+(b0+ay1)*dy1)+jj2*(ax2*dx2+(b0+ay2)*dy2)
enddo

 !Normalise integrals:
circh=f12*circh
circm=f12*circm
dchdt=f12*dchdt

return
end subroutine

!==========================================================================

subroutine cubic(ibeg,iend,dx,dy,b,c)
! Finds the cubic interpolation coefficients for the range of contours
! having nodes ibeg <= i <= iend

implicit none

 !Passed variables:
double precision:: dx(ibeg:iend),dy(ibeg:iend),b(ibeg:iend),c(ibeg:iend)
integer:: ibeg,iend

 !Local variables:
double precision:: u(ibeg:iend),v(ibeg:iend)
double precision:: xx,yy,ds
integer:: i,ia,ib

!------------------------------------------------------------------
 !Compute the contour increments:
do i=ibeg,iend
  ia=nextq(i)
  xx=xq(ia)-xq(i)
  dx(i)=xx-ellx*dble(int(xx*hlxi))
  yy=yq(ia)-yq(i)
  dy(i)=yy-elly*dble(int(yy*hlyi))
enddo

 !Calculate the cubic interpolation coefficients:
do i=ibeg,iend
  v(i)=dx(i)*dx(i)+dy(i)*dy(i)
enddo

do ib=ibeg,iend
  i=nextq(ib)
  u(i)=v(ib)
  b(i)=-dx(ib)
  c(i)=-dy(ib)
enddo

do i=ibeg,iend
  if (dx(i)*b(i)+dy(i)*c(i) .gt. zero) then 
     !Set curvature to zero at corners:
    b(i)=zero
  else
    b(i)=(dx(i)*c(i)-b(i)*dy(i))/ &
     & sqrt((b(i)*v(i)-dx(i)*u(i))**2+(c(i)*v(i)-dy(i)*u(i))**2+small3)
  endif
enddo

do i=ibeg,iend
  ia=nextq(i)
  ds=sqrt(v(i))
  u(i)=ds*(b(ia)+b(i))
  c(i)=ds*(b(ia)-b(i))
enddo

do i=ibeg,iend
  b(i)=f12*(u(i)-c(i))
  c(i)=f13*c(i)
enddo

return
end subroutine

!=======================================================================

 !Main end module
end module
