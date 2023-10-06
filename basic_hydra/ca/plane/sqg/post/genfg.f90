program genfg

! Program to create an image of the buoyancy field on the ultra-fine grid

use contours
use congen

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Declarations:
real:: tr4,qdr4(ny,nx)
character(len=3):: pind
double precision:: qq(ny,nx)

!-----------------------------------------------------------------
write(*,*) ' Initialising ...'
call init_contours

write(*,*) ' What is the numbered suffix of the file to process?'
read(*,*) iind
write(pind(1:3),'(i3.3)') iind

 !Open and read residual needed to build ultra-fine-grid buoyancy:
open(40,file='cont/qqresi.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
read(40,rec=iind) tr4,qdr4
close(40) 
do ix=1,nx
   do iy=1,ny
    qq(iy,ix)=dble(qdr4(iy,ix))
  enddo
enddo

 !Open and read buoyancy contours:
open(40,file='cont/qqsynopsis.asc',status='old')
do i=1,iind
  read(40,*) nq,nptq,tr8,qjump
enddo
close(40)
write(*,*) ' nq = ',nq,' nptq = ',nptq,' t = ',tr8

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

 !Generate and save fine grid field:
call ugridsave(qq,xq,yq,qjump,nextq,indq,npq,i1q,nq,nptq,tr4)

contains 

!======================================================================

subroutine ugridsave(qq,xq,yq,dq,nextq,indq,npq,i1q,nq,nptq,tr4)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(ny,nx)
double precision:: xq(npm),yq(npm)
real tr4
integer:: nextq(npm),indq(nm),npq(nm),i1q(nm)

 !Local arrays:
real:: rqa(nyu,nxu),rrqa(nyu/2,nxu/2)

!----------------------------------------------------------------------
! Perform contour->grid conversion:
call con2ugrid

! Bi-linear interpolate the residual qq to the fine grid and add to qa:
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

 !Convert to real*4:
do ix=1,nxu
  do iy=1,nyu
    rqa(iy,ix)=real(qa(iy,ix))
  enddo
enddo

 !Open output file:
open(44,file='fine/qq'//pind//'.r4',form='unformatted', &
      access='stream',status='replace')

write(*,*)
write(*,*) ' Choose the type of image:'
write(*,*) '   (1) Full domain view using all data points'
write(*,*) '   (2) As above, but using every nth point'
write(*,*) '   (3) Partial domain view'
read(*,*) iopt

write(*,*) ' Writing fine/qq'//pind//'.r4'
if (iopt .eq. 1) then
  nxv=nxu
  nyv=nyu
  write(44) tr4,rqa
else if (iopt .eq. 2) then
  write(*,*) ' Enter the stride, n:'
  read(*,*) inc
  nxv=nxu/inc
  nyv=nyu/inc
  do ix=1,nxv
    ixo=(ix-1)*inc+1
    do iy=1,nyv
      iyo=(iy-1)*inc+1
      rrqa(iy,ix)=rqa(iyo,ixo)
    enddo
  enddo
  write(44) tr4,rrqa(1:nyv,1:nxv)
else
  write(*,'(a,i5,a)') ' There are ',nxu,' x grid points.'
  write(*,*) ' Range of x grid points (i1,i2) on the inversion grid?'
  write(*,*) ' (Note: periodic wrapping is allowed using i2 < i1)'
  read(*,*) ix1,ix2
  ix1v=(ix1-1)*mgu+1
  ix2v=ix2*mgu
   !Allow periodic wrapping:
  nxv=mod(ix2v-ix1v+nxu,nxu)+1

  write(*,'(a,i5,a)') ' There are ',nyu,' y grid points.'
  write(*,*) ' Range of y grid points (i1,i2) on the inversion grid?'
  write(*,*) ' (Note: periodic wrapping is allowed using i2 < i1)'
  read(*,*) iy1,iy2
  iy1v=(iy1-1)*mgu+1
  iy2v=iy2*mgu
   !Allow periodic wrapping:
  nyv=mod(iy2v-iy1v+nyu,nyu)+1
  if (iy1v .lt. iy2v) then
    if (ix1v .lt. ix2v) then
      write(44) tr4,rqa(iy1v:iy2v,ix1v:ix2v)
    else
      write(44) tr4,rqa(iy1v:iy2v,ix1v:nxu),rqa(iy1v:iy2v,1:ix2v)
    endif
  else
    if (ix1v .lt. ix2v) then
      write(44) tr4,(rqa(iy1v:nyu,ix),rqa(1:iy2v,ix),ix=ix1v,ix2v)
    else
      write(44) tr4,(rqa(iy1v:nyu,ix),rqa(1:iy2v,ix),ix=ix1v,nxu), &
                    (rqa(iy1v:nyu,ix),rqa(1:iy2v,ix),ix=1,ix2v)
    endif
  endif
endif

close(44)

write(*,*)
write(*,*) ' Image by typing the command'
write(*,*)
write(*,'(a,i5,1x,i5)') ' dataview fine/qq'//pind//'.r4 -ndim ',nxv,nyv
write(*,*)

return
end subroutine

!==============================================

end program
