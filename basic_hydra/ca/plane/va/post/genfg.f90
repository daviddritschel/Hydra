program genfg

! Program to create an image of the PV field on the ultra-fine grid

use common

implicit none

 !Declarations:
double precision:: q0(ng,ng)
double precision:: t0,t1,delt
integer:: iind,i,j,ibeg,iend
character(len=3):: pind

!-----------------------------------------------------------------
call init_contours

open(40,file='contours/qqsynopsis.asc',status='old')
read(40,*) nq,nptq,t0,qjump
read(40,*) nq,nptq,t1,qjump
rewind(40)
delt=t1-t0

write(*,*)
write(*,*) ' At what time do you wish to create the PV on an ultra-fine grid?'
read(*,*) t
iind=nint(t/delt)+1

do i=1,iind
  read(40,*) nq,nptq,t,qjump
enddo
close(40)
write(*,*)
write(*,'(a,i7,a,i8)') ' number of contours = ',nq, &
                      ', total number of nodes = ',nptq

 !Open and read residual needed to build ultra-fine-grid PV:
open(40,file='contours/qqresi.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(40,rec=iind) t,q0
close(40)

write(pind,'(i3.3)') iind
open(40,file='contours/qqindex'//pind,form='unformatted',status='old')
read(40) npq(1:nq),i1q(1:nq),indq(1:nq)
close(40)

open(40,file='contours/qqnodes'//pind,form='unformatted',status='old')
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
call ugridsave(q0)

contains 

!======================================================================

subroutine ugridsave(q0)

use congen

implicit none

 !Passed array:
double precision:: q0(ng,ng)

 !Local variables:
real:: rqa(ngu,ngu),rrqa(ngu/2,ngu/2)
integer:: ix,ixf,ix0,ix1,ix2,iy,iyf,iy0,iy1,iy2
integer:: iopt,inc,ngv,ixo,iyo,ix1v,iy1v,ix2v,iy2v

!----------------------------------------------------------------------
! Perform contour->grid conversion:
call con2ugrid

! Bi-linear interpolate the residual q0 to the fine grid and add to qa:
do ix=1,ngu
  ixf=ixfw(ix)
  ix0=ix0w(ix)
  ix1=ix1w(ix)

  do iy=1,ngu
    iyf=ixfw(iy)
    iy0=ix0w(iy)
    iy1=ix1w(iy)

    qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*q0(iy0,ix0) &
                     & +w10(iyf,ixf)*q0(iy1,ix0) &
                     & +w01(iyf,ixf)*q0(iy0,ix1) &
                     & +w11(iyf,ixf)*q0(iy1,ix1)
  enddo
enddo

 !Convert to real*4:
rqa=real(qa(1:ngu,1:ngu))

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
  ngv=ngu
  write(44) real(t),rqa
else if (iopt .eq. 2) then
  write(*,*) ' Enter the stride, n:'
  read(*,*) inc
  ngv=ngu/inc
  do ix=1,ngv
    ixo=(ix-1)*inc+1
    do iy=1,ngv
      iyo=(iy-1)*inc+1
      rrqa(iy,ix)=rqa(iyo,ixo)
    enddo
  enddo
  write(44) real(t),rrqa(1:ngv,1:ngv)
else
  write(*,'(a,i5,a)') ' There are ',ngu,' x grid points.'
  write(*,*) ' Range of x grid points (i1,i2) on the inversion grid?'
  write(*,*) ' (Note: periodic wrapping is allowed using i2 < i1)'
  read(*,*) ix1,ix2
  ix1v=(ix1-1)*mgu+1
  ix2v=ix2*mgu
   !Allow periodic wrapping:
  ngv=mod(ix2v-ix1v+ngu,ngu)+1

  write(*,'(a,i5,a)') ' There are ',ngu,' y grid points.'
  write(*,*) ' Range of y grid points (i1,i2) on the inversion grid?'
  write(*,*) ' (Note: periodic wrapping is allowed using i2 < i1)'
  read(*,*) iy1,iy2
  iy1v=(iy1-1)*mgu+1
  iy2v=iy2*mgu
   !Allow periodic wrapping:
  ngv=mod(iy2v-iy1v+ngu,ngu)+1
  if (iy1v .lt. iy2v) then
    if (ix1v .lt. ix2v) then
      write(44) real(t),rqa(iy1v:iy2v,ix1v:ix2v)
    else
      write(44) real(t),rqa(iy1v:iy2v,ix1v:ngu),rqa(iy1v:iy2v,1:ix2v)
    endif
  else
    if (ix1v .lt. ix2v) then
      write(44) real(t),(rqa(iy1v:ngu,ix),rqa(1:iy2v,ix),ix=ix1v,ix2v)
    else
      write(44) real(t),(rqa(iy1v:ngu,ix),rqa(1:iy2v,ix),ix=ix1v,ngu), &
                        (rqa(iy1v:ngu,ix),rqa(1:iy2v,ix),ix=1,ix2v)
    endif
  endif
endif

close(44)

write(*,*)
write(*,*) ' Image by typing the command'
write(*,*)
write(*,'(a,i5,1x,i5,a)') ' dataview fine/qq'//pind//'.r4 -ndim ',ngv,ngv,' &'
write(*,*)

return
end subroutine

!==============================================

end program
