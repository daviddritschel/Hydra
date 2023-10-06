      program vorprop
      implicit double precision (a-h,o-z)

c-----------------------------------------------------------------------
c       Computes vortex properties (area, centroid, deformation) from
c       data in cont/znnn where nnn = 000, 001, ... (selected below) 
c       after running vorid.F

c        *** Adapted by D. G. Dritschel, 19 October 2012 ***

c       Copy into data directory and compile using fc vorprop
c-----------------------------------------------------------------------

c      Specify ng as the grid resolution here:
      parameter (ng=1024,nt=2*ng)
      parameter (ngridp=ng*nt,npm=125*ngridp)
      parameter (nm=npm/20+npm/200)
c      npm:  maximum number of nodes on all contours
c      nm:   maximum number of contours

      parameter (zero=0.d0,one=1.d0,two=2.d0,f12=one/two)
      parameter (f13=one/3.d0,f16=one/6.d0)
      parameter (pi=3.141592653589793238d0,twopi=two*pi,fourpi=4.d0*pi)
      parameter (small=1.d-12,small3=small*small*small)

      parameter (p1=0.21132486541d0,p2=0.78867513459d0)
c      p1,p2: points for 2-point Gaussian Quadrature
      parameter (dp1=two*p1,dp2=two*p2,tp1=1.5d0*p1,tp2=1.5d0*p2)

      dimension x(npm),y(npm),z(npm)
      dimension a(npm),b(npm),c(npm)
      dimension dx(npm),dy(npm),dz(npm)
      dimension ax(npm),ay(npm),az(npm)
      dimension sx(npm),sy(npm),sz(npm)

      dimension xbar(nm),ybar(nm),zbar(nm)
      dimension xvc1(nm),yvc1(nm)
      dimension xvc2(nm),yvc2(nm),zvc2(nm)
      dimension area(nm),def1(nm),def2(nm)

      dimension next(npm),i1(nm),i2(nm),np(nm),ind(nm)
      dimension ivor(nm,2),nvor(2),jvl1(2),jvl2(2)

      character*3 pind,clev

c----------------------------------------------------------
c      Define maximum grid spacing and minimum vortex area:
      dl=twopi/dble(nt)
      aremin=(dl/4.d0)**2

c---------------------------------------------------------------------
c      Select period to analyse:
      write(*,*) ' Which period do you wish to analyse (0 for first):'
      read(*,*) jper
      pind='000'
      write(pind(1:3),'(i3.3)') jper
      loop=jper+1

      write(*,*) ' Select a vorticity level to define vortex edges:'
      read(*,*) indv
      clev='000'
      write(clev(1:3),'(i3.3)') indv

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

c      Define next(i), the node following i:
      do i=1,npt-1
        next(i)=i+1
      enddo
      do j=1,n
        next(i2(j))=i1(j)
      enddo

c      Ensure points all lie on the surface of the unit sphere:
      do i=1,npt
        fac=one/sqrt(x(i)**2+y(i)**2+z(i)**2)
        x(i)=fac*x(i)
        y(i)=fac*y(i)
        z(i)=fac*z(i)
      enddo

c-----------------------------------------------------------
c      Compute cubic coefficients for quadrature below:
      do i=1,npt
        ia=next(i)
        dx(i)=x(ia)-x(i)
        dy(i)=y(ia)-y(i)
        dz(i)=z(ia)-z(i)
        sy(i)=dx(i)**2+dy(i)**2+dz(i)**2
        xtmp=y(i)*dz(i)-z(i)*dy(i)
        ytmp=z(i)*dx(i)-x(i)*dz(i)
        ztmp=x(i)*dy(i)-y(i)*dx(i)
        afac=sqrt(sy(i)/(xtmp**2+ytmp**2+ztmp**2+small))
        ax(i)=afac*xtmp
        ay(i)=afac*ytmp
        az(i)=afac*ztmp
      enddo

      do ib=1,npt
        i=next(ib)
        sx(i)=sy(ib)
        a(i)=-dx(ib)
        b(i)=-dy(ib)
        c(i)=-dz(ib)
      enddo

      do i=1,npt
        if (dx(i)*a(i)+dy(i)*b(i)+dz(i)*c(i) .lt. zero) then
          b(i)=(x(i)*(dy(i)*c(i)-b(i)*dz(i))+
     .          y(i)*(dz(i)*a(i)-c(i)*dx(i))+
     .          z(i)*(dx(i)*b(i)-a(i)*dy(i)))/
     .           sqrt((a(i)*sy(i)-dx(i)*sx(i))**2+
     .                (b(i)*sy(i)-dy(i)*sx(i))**2+
     .                (c(i)*sy(i)-dz(i)*sx(i))**2+small3)
        else
c          Set curvature to zero at corners:
          b(i)=zero
        endif
      enddo

      do i=1,npt
        ia=next(i)
        darc=sqrt(sy(i))
        cadd=darc*(b(ia)+b(i))
        cdif=darc*(b(ia)-b(i))
        a(i)=f16*cdif-f12*cadd
        b(i)=f12*(cadd-cdif)
        c(i)=f13*cdif
      enddo

c-----------------------------------------------------------
c      Read negative vortices:
      m=1
      open(21,file='vorm'//pind//'-'//clev,status='old')
      rewind 21
      read(21,*) nvor(m),indv
      if (nvor(m) .gt. 0) then
        read(21,*) jvl1(m),jvl2(m)
        do j=jvl1(m),jvl2(m)
          read(21,*) ivor(j,m)
        enddo
      endif
      close(21)

c      Read positive vortices:
      m=2
      open(21,file='vorp'//pind//'-'//clev,status='old')
      rewind 21
      read(21,*) nvor(m),indv
      if (nvor(m) .gt. 0) then
        read(21,*) jvl1(m),jvl2(m)
        do j=jvl1(m),jvl2(m)
          read(21,*) ivor(j,m)
        enddo
      endif
      close(21)

      write(*,'(a,f12.5)') ' Analysing data at time t = ',t

      open(40,file='prop'//pind//'-'//clev,status='unknown')
      atot=small
      acoh=zero

c--------------------------------------------------------------
      do m=1,2
        if (nvor(m) .gt. 0) then
c        Compute centroids:
        do k=1,nvor(m)
          xbar(k)=zero
          ybar(k)=zero
          zbar(k)=zero
          xvc2(k)=zero
          yvc2(k)=zero
          zvc2(k)=zero
        enddo

        do j=jvl1(m),jvl2(m)
          k=ivor(j,m)
          do i=i1(j),i2(j)
            sx(i)=dy(i)*az(i)-dz(i)*ay(i)
            sy(i)=dz(i)*ax(i)-dx(i)*az(i)
            sz(i)=dx(i)*ay(i)-dy(i)*ax(i)

c            1st Gaussian point:
            eta=p1*(a(i)+p1*(b(i)+p1*c(i)))
            deta=a(i)+dp1*(b(i)+tp1*c(i))
            del=f12*p1*(one-p1)
            ddel=f12-p1
            x1=x(i)+p1*dx(i)+eta*ax(i)+del*sx(i)
            y1=y(i)+p1*dy(i)+eta*ay(i)+del*sy(i)
            z1=z(i)+p1*dz(i)+eta*az(i)+del*sz(i)
            dx1=dx(i)+deta*ax(i)+ddel*sx(i)
            dy1=dy(i)+deta*ay(i)+ddel*sy(i)
            dz1=dz(i)+deta*az(i)+ddel*sz(i)
            ds1=sqrt(dx1**2+dy1**2+dz1**2)

c            2nd Gaussian point:
            eta=p2*(a(i)+p2*(b(i)+p2*c(i)))
            deta=a(i)+dp2*(b(i)+tp2*c(i))
            del=f12*p2*(one-p2)
            ddel=f12-p2
            x2=x(i)+p2*dx(i)+eta*ax(i)+del*sx(i)
            y2=y(i)+p2*dy(i)+eta*ay(i)+del*sy(i)
            z2=z(i)+p2*dz(i)+eta*az(i)+del*sz(i)
            dx2=dx(i)+deta*ax(i)+ddel*sx(i)
            dy2=dy(i)+deta*ay(i)+ddel*sy(i)
            dz2=dz(i)+deta*az(i)+ddel*sz(i)
            ds2=sqrt(dx2**2+dy2**2+dz2**2)

c            Accumulate integrals:
            xbar(k)=xbar(k)+y1*dz1-z1*dy1+y2*dz2-z2*dy2
            ybar(k)=ybar(k)+z1*dx1-x1*dz1+z2*dx2-x2*dz2
            zbar(k)=zbar(k)+x1*dy1-y1*dx1+x2*dy2-y2*dx2
            xvc2(k)=xvc2(k)+x1*ds1+x2*ds2
            yvc2(k)=yvc2(k)+y1*ds1+y2*ds2
            zvc2(k)=zvc2(k)+z1*ds1+z2*ds2
          enddo
        enddo

c        Normalise and possibly switch direction of centre:
        do k=1,nvor(m)
          ss=sign(one,xbar(k)*xvc2(k)+ybar(k)*yvc2(k)+zbar(k)*zvc2(k))
          fac=ss/sqrt(xbar(k)**2+ybar(k)**2+zbar(k)**2)
          xbar(k)=xbar(k)*fac
          ybar(k)=ybar(k)*fac
          zbar(k)=zbar(k)*fac
c          Define two other unit vectors:
          rbar=sqrt(xbar(k)**2+ybar(k)**2)
          cphi=xbar(k)/rbar
          sphi=ybar(k)/rbar
          xvc1(k)=-sphi
          yvc1(k)= cphi
          xvc2(k)=-cphi*zbar(k)
          yvc2(k)=-sphi*zbar(k)
          zvc2(k)= rbar
        enddo

c        Compute area and deformation integrals:
        do k=1,nvor(m)
          area(k)=zero
          def1(k)=zero
          def2(k)=zero
        enddo

        do j=jvl1(m),jvl2(m)
          k=ivor(j,m)
          do i=i1(j),i2(j)
c            1st Gaussian point:
            eta=p1*(a(i)+p1*(b(i)+p1*c(i)))
            deta=a(i)+dp1*(b(i)+tp1*c(i))
            del=f12*p1*(one-p1)
            ddel=f12-p1
            x1=x(i)+p1*dx(i)+eta*ax(i)+del*sx(i)
            y1=y(i)+p1*dy(i)+eta*ay(i)+del*sy(i)
            z1=z(i)+p1*dz(i)+eta*az(i)+del*sz(i)
            dx1=dx(i)+deta*ax(i)+ddel*sx(i)
            dy1=dy(i)+deta*ay(i)+ddel*sy(i)
            dz1=dz(i)+deta*az(i)+ddel*sz(i)
c            place in coordinates relative to centre of vortex k:
            u1=x1*xvc1(k)+y1*yvc1(k)
            v1=x1*xvc2(k)+y1*yvc2(k)+z1*zvc2(k)
            w1=x1*xbar(k)+y1*ybar(k)+z1*zbar(k)
            du1=dx1*xvc1(k)+dy1*yvc1(k)
            dv1=dx1*xvc2(k)+dy1*yvc2(k)+dz1*zvc2(k)
            da1=(u1*dv1-v1*du1)/(one+w1)
            wc1=one-w1

c            2nd Gaussian point:
            eta=p2*(a(i)+p2*(b(i)+p2*c(i)))
            deta=a(i)+dp2*(b(i)+tp2*c(i))
            del=f12*p2*(one-p2)
            ddel=f12-p2
            x2=x(i)+p2*dx(i)+eta*ax(i)+del*sx(i)
            y2=y(i)+p2*dy(i)+eta*ay(i)+del*sy(i)
            z2=z(i)+p2*dz(i)+eta*az(i)+del*sz(i)
            dx2=dx(i)+deta*ax(i)+ddel*sx(i)
            dy2=dy(i)+deta*ay(i)+ddel*sy(i)
            dz2=dz(i)+deta*az(i)+ddel*sz(i)
            ds2=sqrt(dx2**2+dy2**2+dz2**2)
c            place in coordinates relative to centre of vortex k:
            u2=x2*xvc1(k)+y2*yvc1(k)
            v2=x2*xvc2(k)+y2*yvc2(k)+z2*zvc2(k)
            w2=x2*xbar(k)+y2*ybar(k)+z2*zbar(k)
            du2=dx2*xvc1(k)+dy2*yvc1(k)
            dv2=dx2*xvc2(k)+dy2*yvc2(k)+dz2*zvc2(k)
            da2=(u2*dv2-v2*du2)/(one+w2)
            wc2=one-w2

c            Accumulate integrals:
            area(k)=area(k)+da1+da2
            def1(k)=def1(k)+da1*wc1+da2*wc2
            def2(k)=def2(k)+da1*wc1**2+da2*wc2**2
          enddo
        enddo

c        Normalise:
        do k=1,nvor(m)
          area(k)=abs(f12*area(k))
          def1(k)=abs(def1(k)*area(k)/(def2(k)*twopi))
        enddo
c        Note: def1 = 4*a*(1+a^2)/(3+2*a^2+3*a^4) for a small 
c              elliptical vortex of aspect ratio a.

c        Write data (exclude very small vortices):
        do k=1,nvor(m)
          atot=atot+area(k)
          if (area(k) .gt. aremin) then
            acoh=acoh+area(k)
            write(40,'(i7,5(1x,f13.10))') k,area(k),def1(k),
     .                              xbar(k),ybar(k),zbar(k)
          endif
        enddo
        endif
      enddo

      write(*,'(a,f12.10)') '     Total fractional area = ',atot/fourpi
      write(*,'(a,f12.10)') '  Coherent fractional area = ',acoh/fourpi
      write(*,'(a,f12.10)') '   Coherent/total fraction = ',acoh/atot
      write(*,*)
      write(*,*) 
     .  '  The vortex properties are ready in prop'//pind//'-'//clev

      end
