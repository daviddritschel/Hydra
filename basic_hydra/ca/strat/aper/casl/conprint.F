      program print
c========================================================================
c      Makes a postscript print of contours in the usual CASL format.
c            *** This code is for a non-periodic domain ***

c              -->> need to compile with pfc.
c========================================================================

      parameter (npm=20000000,nm=npm/20)
c      npm,nm: max # of pts and PV contours in data
      parameter (wpage=7.68,nppin=400)
c      wpage: the width (in inches) of the plotted page that 
c             is assumed by the graphics package.
c      nppin: # of pixels per inch on the postscript printer
      parameter (fdash=0.5)
c      fdash: the fractional length of the segment connecting
c             adjacent nodes which is plotted for contours associated
c             with negative PV (having ind(j) < 0).
      parameter (adash=0.5+0.5*fdash,bdash=1.-adash)

      dimension x(npm),y(npm)
      dimension np(nm),i1(nm),ind(nm)
      logical closed(nm)
      character*40 confile

c-----------------------------------------------------------
      write(*,*) ' Enter the name of the contour file:'
      read(*,*) confile
      open(11,file=confile,status='old',err=999)
      rewind(11)

c      Read file containing domain limits:
      open(80,file='domain_limits',status='old')
      read(80,*) xmin,xmax
      read(80,*) ymin,ymax
      close(80)

      write(*,*) '  Width, in inches, of the plot? '
      read(*,*) width

      write(*,*) '  Print just open contours (0) or all contours (1)?'
      read(*,*) iopt

      if (iopt .eq. 0) then
        write(*,*) '  Plot limits, xmin, xmax, ymin, ymax?'
        read(*,*) xminp,xmaxp,yminp,ymaxp
      else
        xminp=xmin
        yminp=ymin
        xmaxp=xmax
        ymaxp=ymax
      endif

      xcen=(xminp+xmaxp)/2.
      ycen=(yminp+ymaxp)/2.
      xhwi=(xmaxp-xminp)/2.
      yhwi=(ymaxp-yminp)/2.

      small=1.e-6*xhwi

      asp=yhwi/xhwi
      height=asp*width
      pscale=xhwi*wpage/width

c      initialize unix graphics:
      call gsetup
      call gscale(-pscale,pscale,-pscale,pscale)
      call izoom
      call arzoom(1.,1.)

      read(11,*) n,npt,t
      do j=1,n
        read(11,*) np(j),i1(j),ind(j)
      enddo
      do i=1,npt
        read(11,*) x(i),y(i)
      enddo

c      Determine closed contours:
      do j=1,n
        ib=i1(j)
        ie=i1(j)+np(j)-1
        closed(j)=((min(abs(x(ib)-xmin),abs(x(ib)-xmax),
     .                  abs(y(ib)-ymin),abs(y(ib)-ymax)) .gt. small)
     .        .or. (min(abs(x(ie)-xmin),abs(x(ie)-xmax),
     .                  abs(y(ie)-ymin),abs(y(ie)-ymax)) .gt. small))
      enddo

c      rotate domain by -90 degrees to make plot come out with
c      its top at the top of the printed page:
      do i=1,npt
        xtmp=x(i)-xcen
        ytmp=y(i)-ycen
        x(i)= ytmp
        y(i)=-xtmp
      enddo

      call gscale(-pscale,pscale,-pscale,pscale)
c      draw a bold border around the frame:
      call ggline(-yhwi,-xhwi, yhwi,-xhwi,3)
      call ggline( yhwi,-xhwi, yhwi, xhwi,3)
      call ggline( yhwi, xhwi,-yhwi, xhwi,3)
      call ggline(-yhwi, xhwi,-yhwi,-xhwi,3)

c      finally, draw the contours:
      if (iopt .eq. 1) then
        do j=1,n
          ib=i1(j)
          ie=i1(j)+np(j)-1
          p=adash+sign(bdash,float(ind(j)))
          do i=ib,ie-1
            call ggline(x(i),y(i),x(i)+p*(x(i+1)-x(i)),
     .                            y(i)+p*(y(i+1)-y(i)),1)
          enddo
c          Connect final and initial point if contour is closed:
          if (closed(j))
     .      call ggline(x(ie),y(ie),x(ie)+p*(x(ib)-x(ie)),
     .                              y(ie)+p*(y(ib)-y(ie)),1)
        enddo
      else
        do j=1,n
          if (.not. closed(j)) then
          ib=i1(j)
          ie=i1(j)+np(j)-1
          p=adash+sign(bdash,float(ind(j)))
          do i=ib,ie-1
            call ggline(x(i),y(i),x(i)+p*(x(i+1)-x(i)),
     .                            y(i)+p*(y(i+1)-y(i)),1)
          enddo
c          Connect final and initial point if contour is closed:
          if (closed(j))
     .      call ggline(x(ie),y(ie),x(ie)+p*(x(ib)-x(ie)),
     .                              y(ie)+p*(y(ib)-y(ie)),1)
          endif
        enddo
      endif

c      terminate graphics:
      call gclear
      call gend
      close(11)
      write(*,*) ' The plot is ready in output.ps'

      stop
999   write(*,*) ' File does not exist!'

      end
