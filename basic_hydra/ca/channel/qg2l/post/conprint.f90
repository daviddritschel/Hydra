program conprint

! Program to output a snapshot of contours extracted form the saved contour file 
! of an aper run to an svg image. Although the image is scalable it has a default
! resolution (xpx,ypx) defined below - note this is used in further conversions to 
! for eg. a png file and should be regarded as the resolution of the image for most 
! purposes.

use contours
use congen

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Declarations:
integer,parameter:: xpx=nxu,ypx=nyu
integer,parameter:: nskip=1
character(len=3):: pind
integer:: iclosed(nm)

!-----------------------------------------------------------------
write(*,*) ' Preparing...'
call init_contours

write(*,*) ' Which layer do you want to process?'
read(*,*) layt

write(*,*) ' What is the numbered suffix of the file?'
read(*,*) iind
write(pind(1:3),'(i3.3)') iind

 !Open contour files:
open(40,file='contours/qqsynopsis.asc',status='old')
do i=1,iind
  read(40,*) nq,nptq,t,qjump,qavg
enddo
close(40)

open(40,file='contours/qqindex'//pind,form='unformatted', &
    & access='direct',status='old',recl=20*nq)
read(40,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq),layq(1:nq),iclosed(1:nq)
close(40)

open(40,file='contours/qqnodes'//pind,form='unformatted', &
    & access='direct',status='old',recl=16*nptq)
read(40,rec=1) xq(1:nptq),yq(1:nptq)
close(40)

 !Work out ranges of contours in chosen layer:
do iz=1,nz
  jl2q(iz)=0
enddo

do j=1,nq
  i2q(j)=i1q(j)+npq(j)-1
  il2q(layq(j))=i2q(j)
  jl2q(layq(j))=j
enddo

jbeg=1
do iz=1,nz
  if (jl2q(iz) .gt. 0) then
    jl1q(iz)=jbeg
    jbeg=jl2q(iz)+1
  endif
enddo

 !Reconstruct nextq array:
do j=jl1q(layt),jl2q(layt)
  ibeg=i1q(j)
  iend=ibeg+npq(j)-1
  do i=ibeg,iend-1
    nextq(i)=i+1
  enddo
  if (iclosed(j) .eq. 0) then
    nextq(iend)=0
  else
    nextq(iend)=ibeg
  endif
enddo 

call conwrite(layt)

contains 

!-----------------------------------------------------------------
subroutine conwrite(iz)

! Write an svg file directly with a path tracing each individual contour
! in layer iz

implicit double precision(a-h,o-z)
implicit integer(i-n)

!-------------------------------------------------------------------------
open(88,file='fine/out.svg',status='unknown')

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
do j=jl1q(iz),jl2q(iz),nskip
  i1=i1q(j)
  write(88,'(a,2(f8.1))') '<path d="M',(xq(i1)-xmin)*scalex,(ymax-yq(i1))*scaley 
  i=i1
  ia=nextq(i1)
  do while (ia .ne. 0)
    if (abs(xq(ia)-xq(i)) .lt. hlx) then 
      i=ia
      write(88,'(a,2(f8.1))') 'L',(xq(i)-xmin)*scalex,(ymax-yq(i))*scaley 
    else
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
    endif
    if (i .eq. i1) exit
    ia=nextq(i)   
  enddo
   !Close the contour write:
  write(88,'(a)') ' "      fill="none" stroke="black" stroke-width="0.5"/>'
enddo

 !Write footer for output file:
!---------------------include and edit for a box around the image--------------------
!write(88,'(a,i5,a,i5,a)') '  <rect x="0" y="0" width="',xpx-1,'" height="',ypx-1,'"'
!write(88,'(a)') '        fill="none" stroke="black" stroke-width="1.0" />'
!------------------------------------------------------------------------------------
write(88,'(a)') '</svg>'
close(88)

return
end subroutine

!==============================================

end program
