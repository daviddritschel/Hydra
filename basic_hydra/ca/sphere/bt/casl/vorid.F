      program vorid
      implicit double precision (a-h,o-z)
c-----------------------------------------------------------------------
c       Vortex Identification Program for caps in spherical geometry

c       Finds 2D contiguous vortices from data in cont/znnn where
c       nnn = 000, 001, ... (selected below)

c        *** Adapted by D. G. Dritschel, 16 October 2012 ***

c       Contours (enclosing less than half of the sphere's surface 
c       area) are assumed to be right-handed unless they are holes
c       within other contours.

c       The output is in two files, vorp.dat & vorm.dat, for 
c       positive and negative vortices respectively.  These files
c       contain a listing of the contours belonging to each vortex 
c       in the following format:

c         nvor, lev <--- the total # of vortices found and the vorticity 
c                        level, and then (if nvor > 0), 
c         jvl1,jvl2 <--- the beginning & ending contours C_j having
c                        ind(j) = +/-indv; NB: jvl2 = 0 means
c                        that there are no contours with this index;
c                        and then (if jvl2 > 0), for j = jvl1, ..., jvl2,
c         ivor      <--- the vortex to which C_j belongs.  If zero,
c                        it means that C_j is an improper contour
c                        (see notes in subroutine getvor below).

c       These data can be used to compute diagnostics such as vortex mean 
c       radii, etc, in the separate programme vorprop.F

c       Copy into data directory and compile using fc vorid
c-----------------------------------------------------------------------
c      Specify ng as the grid resolution here and in the subroutine below:
      parameter (ng=1024,nt=2*ng)
      parameter (ngridp=ng*nt,npm=125*ngridp)
      parameter (nm=npm/20+npm/200)
c      npm:  maximum number of nodes on all contours
c      nm:   maximum number of contours

      parameter (zero=0.d0,one=1.d0,two=2.d0)

      common /nodes/ x(npm),y(npm),z(npm)
      common /cattr/ i1(nm),i2(nm),np(nm),ind(nm)
      common /vortx/ ivor(nm),jvl1,jvl2
      common /count/ nvor,indv,n,npt

      character*3 pind,clev
c---------------------------------------------------------------------
c      Select period to analyse:
      write(*,*) ' Which period do you wish to analyse (0 for first):'
      read(*,*) jper
      pind='000'
      write(pind(1:3),'(i3.3)') jper
      loop=jper+1

c      Read vorticity contours:
      open(40,file='cont/synopsis.asc',status='old')
      do i=1,loop
        read(40,*) t,n,npt
      enddo
      close(40)
c      n:   the total number of contours at time t
c      npt: the total number of nodes at time t

      open(40,file='cont/index'//pind,form='unformatted',status='old')
        read(40) np(1:n),i1(1:n),ind(1:n)
      close(40)
c      np(j): the number of nodes on C_j
c      i1(j): the index of the first node on C_j
c      ind(j): the field isolevel label associated with C_j

      open(40,file='cont/nodes'//pind,form='unformatted',status='old')
        read(40) x(1:npt),y(1:npt),z(1:npt)
      close(40)
c      x(i),y(i),z(i): the coordinates of node i on a contour.
c      For a given contour j, the nodes range from i1(j) to 
c      i1(j)+np(j)-1 ( = i2(j) just below).

      do j=1,n
        i2(j)=i1(j)+np(j)-1
      enddo

c      Ensure points all lie on the surface of the unit sphere:
      do i=1,npt
        fac=one/sqrt(x(i)**2+y(i)**2+z(i)**2)
        x(i)=fac*x(i)
        y(i)=fac*y(i)
        z(i)=fac*z(i)
      enddo

      write(*,'(a,f12.5)') ' Analysing data at time t = ',t
      write(*,*)
      write(*,*) ' Select a vorticity level to define vortex edges:'
      read(*,*) indv
      clev='000'
      write(clev(1:3),'(i3.3)') indv

c--------------------------------------------------------------
c      Get negative field level vortices having ind(j)=-|indv|:
      indv=-indv
      nvor=0
      call getvor

c      Write data:
      open(21,file='vorm'//pind//'-'//clev,status='unknown')
      rewind 21
      write(21,'(i7,1x,i5)') nvor,indv
      if (nvor .gt. 0) then
        write(21,'(i7,1x,i7)') jvl1,jvl2
        if (jvl2 .gt. 0) then
          do j=jvl1,jvl2
            write(21,'(i7)') ivor(j)
          enddo
        endif
      endif
      close(21)

c--------------------------------------------------------------
c      Get positive field level vortices having ind(j)=+|indv|:
      indv=-indv
      nvor=0
      call getvor

c      Write data:
      open(22,file='vorp'//pind//'-'//clev,status='unknown')
      rewind 22
      write(22,'(i7,1x,i5)') nvor,indv
      if (nvor .gt. 0) then
        write(22,'(i7,1x,i7)') jvl1,jvl2
        if (jvl2 .gt. 0) then
          do j=jvl1,jvl2
            write(22,'(i7)') ivor(j)
          enddo
        endif
      endif
      close(22)

      write(*,*) '  The data are ready in vorp'//pind//'-'//clev//
     .   ' & vorm'//pind//'-'//clev
      write(*,*) 'for positive and negative vortices, respectively.'

      end

c=====================================================================

      subroutine getvor
      implicit double precision(a-h,o-z)
c      Finds all the contours C_j belonging to each vortex; all
c      the C_j have ind(j) = indv, i.e. a given field isolevel.

c      Specify ng as the grid resolution here and in the main programme above
      parameter (ng=1024,nt=2*ng)
      parameter (ngridp=ng*nt,npm=125*ngridp)
      parameter (nm=npm/20+npm/200)
c      npm:  maximum number of nodes on all contours
c      nm:   maximum number of contours
      common /nodes/ x(npm),y(npm),z(npm)
      common /cattr/ i1(nm),i2(nm),np(nm),ind(nm)
      common /vortx/ ivor(nm),jvl1,jvl2
      common /count/ nvor,indv,n,npt

c      Local parameters and arrays:
      parameter (zero=0.d0,one=1.d0,two=2.d0,half=one/two)
      parameter (pi=3.141592653589793238d0,twopi=two*pi)
      parameter (small=1.d-13)

      dimension u(npm),v(npm),w(npm)
      dimension xbar(nm),ybar(nm),zbar(nm)
      dimension xvc1(nm),yvc1(nm)
      dimension xvc2(nm),yvc2(nm),zvc2(nm)
      dimension area(nm),arem(nm),cthe(nm),sthe(nm)
      dimension jreg(nm),jpos(nm),kpos(nm)
      logical bound(nm)

c-------------------------------------------------------------------
c      Get beginning and ending contours (jvl1,jvl2) having 
c      selected index value (ind(j) = indv):
      jvl1=0
      jvl2=0
      do j=1,n
        bound(j)=(ind(j) .eq. indv)
      enddo
      do j=1,n
        if (bound(j)) then
          jvl2=j
          if (jvl1 .eq. 0) jvl1=j
        endif
      enddo
c      If there are no contours with ind(j) = indv, there is nothing to do:
      if (jvl2 .eq. 0) return

      write(*,'(a,i3,a,i7,a)') ' Processing level ',indv,' having ',
     .             jvl2-jvl1+1,' contours'

c-----------------------------------------------------------------------
c      For each region R_j, bounded by a contour C_j with positive area,
c      determine which contours, if any, are contained within it.  These
c      interior contours must have negative area and are thus holes in R_j.

c      Calculate contour centres (xbar,ybar,zbar) to determine handedness:
      do j=jvl1,jvl2
        ibeg=i1(j)
        iend=i2(j)
        xs=x(iend)+x(ibeg)
        ys=y(iend)+y(ibeg)
        zs=z(iend)+z(ibeg)
        xd=x(ibeg)-x(iend)
        yd=y(ibeg)-y(iend)
        zd=z(ibeg)-z(iend)
        xbar(j)=ys*zd-zs*yd
        ybar(j)=zs*xd-xs*zd
        zbar(j)=xs*yd-ys*xd
c        Also need a crude estimate to ensure correct hemisphere:
        sd=sqrt(xd**2+yd**2+zd**2)
        xvc2(j)=xs*sd
        yvc2(j)=ys*sd
        zvc2(j)=zs*sd
      enddo
c      The above initialises (xbar,ybar,zbar) with one segment's contribution
      do j=jvl1,jvl2
        do i=i1(j),i2(j)-1
          xs=x(i)+x(i+1)
          ys=y(i)+y(i+1)
          zs=z(i)+z(i+1)
          xd=x(i+1)-x(i)
          yd=y(i+1)-y(i)
          zd=z(i+1)-z(i)
          xbar(j)=xbar(j)+ys*zd-zs*yd
          ybar(j)=ybar(j)+zs*xd-xs*zd
          zbar(j)=zbar(j)+xs*yd-ys*xd
          sd=sqrt(xd**2+yd**2+zd**2)
          xvc2(j)=xvc2(j)+xs*sd
          yvc2(j)=yvc2(j)+ys*sd
          zvc2(j)=zvc2(j)+zs*sd
        enddo
      enddo
c      Construct also two unit vectors perpendicular to (xbar,ybar,zbar):
      do j=jvl1,jvl2
        fac=one/sqrt(xbar(j)**2+ybar(j)**2+zbar(j)**2)
        xbar(j)=xbar(j)*fac
        ybar(j)=ybar(j)*fac
        zbar(j)=zbar(j)*fac
        fac=one/sqrt(xvc2(j)**2+yvc2(j)**2+zvc2(j)**2)
        xvc2(j)=xvc2(j)*fac
        yvc2(j)=yvc2(j)*fac
        zvc2(j)=zvc2(j)*fac
c        Possibly correct centre to place in correct hemisphere:
        ss=sign(one,xbar(j)*xvc2(j)+ybar(j)*yvc2(j)+zbar(j)*zvc2(j))
        xbar(j)=xbar(j)*ss
        ybar(j)=ybar(j)*ss
        zbar(j)=zbar(j)*ss
c        Define two other unit vectors:
        rbar=sqrt(xbar(j)**2+ybar(j)**2)
        cphi=xbar(j)/rbar
        sphi=ybar(j)/rbar
        xvc1(j)=-sphi
        yvc1(j)= cphi
        xvc2(j)=-cphi*zbar(j)
        yvc2(j)=-sphi*zbar(j)
        zvc2(j)= rbar
      enddo
c      Define coordinates of C_j, (u,v,w), relative to centre of C_j:
      do j=jvl1,jvl2
        do i=i1(j),i2(j)
          u(i)=x(i)*xvc1(j)+y(i)*yvc1(j)
          v(i)=x(i)*xvc2(j)+y(i)*yvc2(j)+z(i)*zvc2(j)
          w(i)=x(i)*xbar(j)+y(i)*ybar(j)+z(i)*zbar(j)
        enddo
      enddo

c      Calculate contour areas to determine handedness:
      do j=jvl1,jvl2
        ibeg=i1(j)
        iend=i2(j)
        xs=half*(u(iend)+u(ibeg))
        ys=half*(v(iend)+v(ibeg))
        zs=half*(w(iend)+w(ibeg))
        xd=u(ibeg)-u(iend)
        yd=v(ibeg)-v(iend)
        zd=w(ibeg)-w(iend)
        area(j)=(xs*yd-ys*xd)/(one+zs)
      enddo
c      The above initialises area(j) with one segment's contribution
      do j=jvl1,jvl2
        do i=i1(j),i2(j)-1
          xs=half*(u(i)+u(i+1))
          ys=half*(v(i)+v(i+1))
          zs=half*(w(i)+w(i+1))
          xd=u(i+1)-u(i)
          yd=v(i+1)-v(i)
          zd=w(i+1)-w(i)
          area(j)=area(j)+(xs*yd-ys*xd)/(one+zs)
        enddo
      enddo

c      Include sign of vorticity jump in area(j) (to identify bounding
c      contours:
      ss=sign(one,dble(indv))
      do j=jvl1,jvl2
        area(j)=ss*area(j)
        bound(j)=(area(j) .gt. zero)
c        contours C_j with positive area bound the regions R_j.
        if (bound(j)) then
          jreg(j)=j
        else
          jreg(j)=0
c          jreg gives the region to which C_j is associated; 
c          it is computed immediately below.
        endif
      enddo

c      Next determine which contours are contained within which.
c      To speed up the computation, the bounding circle of each
c      contour is first calculated:
      do j=jvl1,jvl2
        cthe(j)=one
        do i=i1(j),i2(j)
          cthe(j)=min(cthe(j),w(i))
        enddo
        sthe(j)=sqrt(one-cthe(j)**2)
      enddo
c      cthe=cos(theta), sthe=sin(theta): theta is the solid angular
c      radius of the circle centred at (xbar,ybar,zbar) enclosing C_j.

c      Now loop over contour pairs to see which negative area 
c      contours are inside those with positive area:
      npos=0
c      npos below will count the number of cases when C_k has been
c      found to lie within more than one region.
      do j=jvl1,jvl2
        if (bound(j)) then
c          Only look for contours C_k that are inside contours C_j
c          with positive area:
          ibeg=i1(j)
          iend=i2(j)
          do 10 k=jvl1,jvl2
            if (bound(k)) goto 10
c            Only consider negative area contours C_k.  Reduce search
c            to the C_k whose bounding circle overlaps that of C_j:
            if (xbar(j)*xbar(k)+ybar(j)*ybar(k)+zbar(j)*zbar(k) .lt.
     .          cthe(j)*cthe(k)-sthe(j)*sthe(k)) goto 10
            if (k .eq. j) goto 10
c            For C_k to be contained within C_j, require that MOST of the
c            points ik on C_k lie within C_j; check using winding number:
            swind=zero
            do ik=i1(k),i2(k)
c              Define polar stereographic projection of points ij on C_j
c              relative to the point ik on C_k:
              z2=sqrt(x(ik)**2+y(ik)**2)
              y1= x(ik)/z2
              x1=-y(ik)/z2
              x2=-y1*z(ik)
              y2= x1*z(ik)
              do ij=ibeg,iend
                pfac=one/(one+x(ij)*x(ik)+y(ij)*y(ik)+z(ij)*z(ik))
                u(ij)=pfac*(x(ij)*x1+y(ij)*y1)
                v(ij)=pfac*(x(ij)*x2+y(ij)*y2+z(ij)*z2)
              enddo
c              Initialise winding number with first segment's contribution:
              delx=u(iend)
              dely=v(iend)
              xd=u(ibeg)-u(iend)
              yd=v(ibeg)-v(iend)
              wind=atan2(u(iend)*v(ibeg)-v(iend)*u(ibeg),
     .             small+u(iend)*u(ibeg)+v(iend)*v(ibeg))
              do ij=ibeg,iend-1
                delx=u(ij)-up
                dely=v(ij)-vp
                xd=u(ij+1)-u(ij)
                yd=v(ij+1)-v(ij)
                wind=wind+atan2(u(ij)*v(ij+1)-v(ij)*u(ij+1),
     .                    small+u(ij)*u(ij+1)+v(ij)*v(ij+1))
              enddo
              swind=swind+wind
            enddo

c            Examine average winding number to determine topology:
            if (nint(ss*swind/(twopi*dble(np(k)))) .eq. 1) then
c            C_k is (considered to be) inside C_j if the above holds.
              if (jreg(k) .eq. 0) then
                jreg(k)=j
              else
c              C_k has already been found to lie within another
c              region; record k and j for finding out below
c              which is the correct region C_k belongs to:
                npos=npos+1
                kpos(npos)=k
                jpos(npos)=j
              endif
            endif
10        continue
        endif
      enddo
c      Any C_k with jreg(k) = 0 at this stage is an error; e.g.
c      a left-handed contour not contained within another.  This
c      may occur as a result of surgery or node redistribution, 
c      but it is a rare event and normally results in 
c      surgical-scale contours of little significance.

      if (npos .gt. 0) then
c      Process the list of contours that had been assigned to 
c      more than one region:
        do ipos=1,npos
          k=kpos(ipos)
          arem(k)=area(jreg(k))
        enddo
        do ipos=1,npos
          k=kpos(ipos)
          j=jpos(ipos)
          if (area(j) .lt. arem(k)) then
c          the positive-area contour C_j of least area is the parent 
c          of the inner negative-area contour C_k.  C_k is inside R_j.
            arem(k)=area(j)
            jreg(k)=j
          endif
        enddo
      endif

c-------------------------------------------------------------------
c      Identify all contours belonging to a given vortex:

      do j=jvl1,jvl2
        if (bound(j)) then
c        Above, contours with bound(j)=.true. bound the regions R_j.
          nvor=nvor+1
          ivor(j)=nvor
        endif
      enddo

c      Here, we record all contours belonging to a given vortex using
c      the fact that jreg(k) = j means that C_k belongs to R_j:
      do k=jvl1,jvl2
        if (.not. bound(k)) then
          j=jreg(k)
          if (j .gt. 0) then
            ivor(k)=ivor(j)
          else
            ivor(k)=0
            write(*,'(a,i7,a,f12.9)') ' Contour ',k,
     .            ' not associated; area = ',area(k)
          endif
        endif
      enddo
c      The above states that C_k is part of same vortex as its
c      parent region R_j, where j = jreg(k).  Error contours with
c      jreg(k)=0 are assigned a zero vortex index.
      return
      end
