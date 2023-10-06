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
character(len=3):: pind
integer:: iclosed(nm)

!-----------------------------------------------------------------
write(*,*) ' Preparing...'
call init_contours

write(*,*) ' Which set of contours do you want to process'
write(*,*) ' Choose (0) buoyancy or (1) vorticity:'
read(*,*) iopt

write(*,*) ' What is the numbered suffix of the file?'
read(*,*) iind
write(pind(1:3),'(i3.3)') iind


if (iopt .eq. 0) then
   !Open buoyancy contours:
  open(40,file='cont/bbsynopsis.asc',status='old')
  do i=1,iind
    read(40,*) nb,nptb,t,bjump,bavg
  enddo
  close(40)

  open(40,file='cont/bbindex'//pind,form='unformatted', &
      & access='direct',status='old',recl=16*nb)
  read(40,rec=1) npb(1:nb),i1b(1:nb),indb(1:nb),iclosed(1:nb)
  close(40)

  open(40,file='cont/bbnodes'//pind,form='unformatted', &
      & access='direct',status='old',recl=16*nptb)
  read(40,rec=1) xb(1:nptb),yb(1:nptb)
  close(40)

   !Reconstruct nextb array:
  do j=1,nb
    ibeg=i1b(j)
    iend=ibeg+npb(j)-1
    do i=ibeg,iend-1
      nextb(i)=i+1
    enddo
    if (iclosed(j) .eq. 0) then 
      nextb(iend)=0
    else
      nextb(iend)=ibeg
    endif
  enddo 

  call conwrite(xb,yb,bjump,bavg,nextb,indb,npb,i1b,nb,nptb,t,0)

else
   !Open vorticity contours:
  open(40,file='cont/zzsynopsis.asc',status='old')
  do i=1,iind
    read(40,*) nz,nptz,t,zjump,zavg
  enddo
  close(40)

  if (nz .gt. 0) then
    open(40,file='cont/zzindex'//pind,form='unformatted', &
        & access='direct',status='old',recl=16*nz)
    read(40,rec=1) npz(1:nz),i1z(1:nz),indz(1:nz),iclosed(1:nz)
    close(40)

    open(40,file='cont/zznodes'//pind,form='unformatted', &
        & access='direct',status='old',recl=16*nptz)
    read(40,rec=1) xz(1:nptz),yz(1:nptz)
    close(40)

     !Reconstruct nextz array:
    do j=1,nz
      ibeg=i1z(j)
      iend=ibeg+npz(j)-1
      do i=ibeg,iend-1
        nextz(i)=i+1
      enddo
      if (iclosed(j) .eq. 0) then
        nextz(iend)=0
      else
        nextz(iend)=ibeg
      endif
    enddo 
  endif

  call conwrite(xz,yz,zjump,zavg,nextz,indz,npz,i1z,nz,nptz,t,1)
 
endif


contains 

!-----------------------------------------------------------------
subroutine conwrite(xq,yq,dq,qavg,nextq,indq,npq,i1q,nq,nptq,t,ires)

! Write an svg file directly with a path tracing each individual contour.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
integer,parameter:: dbleint=selected_int_kind(16)
integer(kind=dbleint),parameter:: npu=int(nxu+1,kind=dbleint)*int(nyu+1,kind=dbleint)
double precision:: xq(npm),yq(npm)
integer:: nextq(npm),indq(nm),npq(nm),i1q(nm)

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
do j=1,nq
  i1=i1q(j)
  write(88,'(a,2(f8.1))') '<path d="M',(xq(i1)-xmin)*scalex,(ymax-yq(i1))*scaley 
  i=nextq(i1)
  do while ((i .ne. i1) .and. (i .ne. 0))
    write(88,'(a,2(f8.1))') 'L',(xq(i)-xmin)*scalex,(ymax-yq(i))*scaley 
    i=nextq(i)
  enddo 
   !Optionally close the contour:
  if (i .eq. 0 ) then 
    write(88,'(a)') ' "      fill="none" stroke="black" stroke-width="0.5"/>'
  else
    write(88,'(a)') ' z"      fill="none" stroke="black" stroke-width="0.5"/>'
  endif
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
