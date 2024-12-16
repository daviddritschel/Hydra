program conprint

! Program to output a snapshot of contours extracted form the saved contour file 
! of an aper run to an svg image. Although the image is scalable it has a default
! resolution (xpx,ypx) defined below - note this is used in further conversions to 
! for eg. a png file and should be regarded as the resolution of the image for most 
! purposes.

use contours

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Declarations:
integer,parameter:: xpx=nxu,ypx=nyu
integer,parameter:: nskip=1
character(len=3):: pind

!-----------------------------------------------------------------
write(*,*) ' Preparing...'
call init_contours

write(*,*) ' What is the numbered suffix of the file?'
read(*,*) iind
write(pind(1:3),'(i3.3)') iind

 !Open vorticity contours:
open(40,file='cont/qqsynopsis.asc',status='old')
do i=1,iind
  read(40,*) nq,nptq,t,qjump
enddo
close(40)

write(*,*) nq,nptq,t,qjump

if (nq .gt. 0) then
  open(40,file='cont/qqindex'//pind,form='unformatted',status='old')
  read(40) npq(1:nq),i1q(1:nq),indq(1:nq)
  close(40)

  open(40,file='cont/qqnodes'//pind,form='unformatted',status='old')
  read(40) xq(1:nptq),yq(1:nptq)
  close(40)

     !Reconstruct nextq array:
  do j=1,nq
    ibeg=i1q(j)
    iend=ibeg+npq(j)-1
    do i=ibeg,iend-1
      nextq(i)=i+1
    enddo
    nextq(iend)=ibeg
  enddo 
endif

call conwrite

contains 

!-----------------------------------------------------------------
subroutine conwrite

! Write an svg file directly with a path tracing each individual contour.

implicit double precision(a-h,o-z)
implicit integer(i-n)

logical:: crossx,crossy
!-------------------------------------------------------------------------

if (nq .eq. 0) then
  write(*,*) ' No contours read, stopping....'
  stop
endif

open(88,file='out.svg',status='unknown')

 !Set scale factors for fitting contours to the plot domain:
scalex=dble(xpx)/ellx
scaley=dble(ypx)/elly

 !Write header for output file:
write(88,'(a)') '<?xml version="1.0" standalone="no"?>'
write(88,'(a)') '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" '
write(88,'(a)') ' "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">'
write(88,'(a,i5,a,i5,a)') '<svg width="',xpx,'" height="',ypx,'" version="1.1"'
write(88,'(a)') '    xmlns="http://www.w3.org/2000/svg">'

 !For each contour trace out a path:
do j=1,nq,nskip
  i1=i1q(j)
  write(88,'(a,2(f8.1))') '<path d="M',(xq(i1)-xmin)*scalex,(ymax-yq(i1))*scaley 
  i=i1
  ia=nextq(i1)
  do
    crossx=(abs(xq(ia)-xq(i)) .gt. hlx)
    crossy=(abs(yq(ia)-yq(i)) .gt. hly)
    if ((.not. crossx) .and. (.not. crossy)) then 
      i=ia
      write(88,'(a,2(f8.1))') 'L',(xq(i)-xmin)*scalex,(ymax-yq(i))*scaley 
    else 
      if (crossx .and. crossy) then
!        write(*,*) 'Both crossx and crossy true'
        xxa=sign(hlx,xq(i))
        xxb=sign(hlx,xq(ia))
        aa=abs(xxa-xq(i))
        bb=abs(xxb-xq(ia))
        px=aa/(aa+bb)
        yya=sign(hly,yq(i))
        yyb=sign(hly,yq(ia))
        aa=abs(yya-yq(i))
        bb=abs(yyb-yq(ia))
        py=aa/(aa+bb)
        if (px .lt. py) then
!          write(*,*) 'px < py',px,py
          ytm=yq(ia)-yq(i)
          dy=ytm-elly*int(ytm*hlyi)
          yyy=yq(i)+px*dy
          xtm=xq(ia)-xq(i)
          dx=xtm-ellx*int(xtm*hlxi)
          xtm=xq(i)+py*dx
          xxx=xtm-ellx*int(xtm*hlxi)
           !Write out relevant line segments:    
          i=ia
          write(88,'(a,2(f8.1))') 'L',(xxa-xmin)*scalex,(ymax-yyy)*scaley 
          write(88,'(a,2(f8.1))') 'M',(xxb-xmin)*scalex,(ymax-yyy)*scaley 
          write(88,'(a,2(f8.1))') 'L',(xxx-xmin)*scalex,(ymax-yya)*scaley 
          write(88,'(a,2(f8.1))') 'M',(xxx-xmin)*scalex,(ymax-yyb)*scaley 
          write(88,'(a,2(f8.1))') 'L',(xq(i)-xmin)*scalex,(ymax-yq(i))*scaley 
        else
!          write(*,*) 'px > py',px,py
          xtm=xq(ia)-xq(i)
          dx=xtm-ellx*int(xtm*hlxi)
          xxx=xq(i)+py*dx
          ytm=yq(ia)-yq(i)
          dy=ytm-elly*int(ytm*hlyi)
          ytm=yq(i)+px*dy
          yyy=ytm-elly*int(ytm*hlyi)
           !Write out relevant line segments:    
          i=ia
          write(88,'(a,2(f8.1))') 'L',(xxx-xmin)*scalex,(ymax-yya)*scaley 
          write(88,'(a,2(f8.1))') 'M',(xxx-xmin)*scalex,(ymax-yyb)*scaley 
          write(88,'(a,2(f8.1))') 'L',(xxa-xmin)*scalex,(ymax-yyy)*scaley 
          write(88,'(a,2(f8.1))') 'M',(xxb-xmin)*scalex,(ymax-yyy)*scaley 
          write(88,'(a,2(f8.1))') 'L',(xq(i)-xmin)*scalex,(ymax-yq(i))*scaley 
        endif
      else 
        if (crossx) then
!          write(*,*) 'Only crossx true'
           !Wrapping contour - first interpolate endpoints at domain edges:
          xxa=sign(hlx,xq(i))
          xxb=sign(hlx,xq(ia))
          afac=abs(xxa-xq(i))
          bfac=abs(xxb-xq(ia))
          p=afac/(afac+bfac)
          yyy=yq(i)+p*(yq(ia)-yq(i))      
           !Write out relevant line segments:    
          i=ia
          write(88,'(a,2(f8.1))') 'L',(xxa-xmin)*scalex,(ymax-yyy)*scaley 
          write(88,'(a,2(f8.1))') 'M',(xxb-xmin)*scalex,(ymax-yyy)*scaley 
          write(88,'(a,2(f8.1))') 'L',(xq(i)-xmin)*scalex,(ymax-yq(i))*scaley 
        else
!          write(*,*) 'Only crossy true'
           !Wrapping contour - first interpolate endpoints at domain edges:
          yya=sign(hly,yq(i))
          yyb=sign(hly,yq(ia))
          afac=abs(yya-yq(i))
          bfac=abs(yyb-yq(ia))
          p=afac/(afac+bfac)
          xxx=xq(i)+p*(xq(ia)-xq(i))      
           !Write out relevant line segments:    
          i=ia
          write(88,'(a,2(f8.1))') 'L',(xxx-xmin)*scalex,(ymax-yya)*scaley 
          write(88,'(a,2(f8.1))') 'M',(xxx-xmin)*scalex,(ymax-yyb)*scaley 
          write(88,'(a,2(f8.1))') 'L',(xq(i)-xmin)*scalex,(ymax-yq(i))*scaley
        endif
      endif
    endif
    if (i .eq. i1) exit
    ia=nextq(i)   
  enddo
   !Close the contour write:
  write(88,'(a)') ' "      fill="none" stroke="black" stroke-width="0.5"/>'
enddo

 !Write on axes:
!---------------------include if you want axes--------------------
!write(88,'(a,2(f8.1))') '<path d="M',(-ellx-xmin)*scalex,(ymax-zero)*scaley 
!write(88,'(a,2(f8.1))') 'L',(ellx-xmin)*scalex,(ymax-zero)*scaley
!write(88,'(a,2(f8.1))') 'M',(zero-xmin)*scalex,(ymax-elly)*scaley 
!write(88,'(a,2(f8.1))') 'L',(zero-xmin)*scalex,(ymax+elly)*scaley
!write(88,'(a)') ' "      fill="none" stroke="black" stroke-width="0.5"/>'
!------------------------------------------------------------------
 !Write footer for output file:
!--------------include and edit for a box around the image--------------------
!write(88,'(a,i5,a,i5,a)') '  <rect x="0" y="0" width="',xpx-1,'" height="',ypx-1,'"'
!write(88,'(a)') '        fill="none" stroke="black" stroke-width="1.0" />'
!-----------------------------------------------------------------------------
write(88,'(a)') '</svg>'
close(88)

return
end subroutine

!==============================================

end program
