program genfg

! Program for generating an ultra-fine grid image of the PV field in each 
! layer from the contour and residual files created by a casl simulation.

use congen

implicit none

integer:: iind

!-----------------------------------------------------------------
write(*,*) ' Preparing...'
call init_contours

write(*,*) ' What is the numbered suffix of the file?'
read(*,*) iind

 !Read data and initialise contour arrays as well as gridded residual PV:
call readcont(iind,0)

call ugridsave(qd)

contains 

!==========================================================================
subroutine ugridsave(qq)

implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: dbleint=selected_int_kind(16)
integer(kind=dbleint),parameter:: npu=int(nxu,kind=dbleint)*int(nyu+1,kind=dbleint)
real:: qar4(0:nyu,0:nxum1)
double precision:: qq(0:ny,0:nxm1,nz)
character(len=3):: pind

 !------------------------------------------------------------------------
 !Open output file:
write(pind(1:3),'(i3.3)') iind
open(44,file='fine/qq'//pind//'.r4',form='unformatted', &
      access='direct',status='replace',recl=4+npu*4)

do iz=1,nz
   !Convert contours to gridded values (if there are contours in this layer):
  if (jl2q(iz) .gt. 0) call con2ugrid(iz)

   !Bi-linear interpolate the residual qq to the fine grid and add to qa:
  do ix=0,nxum1
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    do iy=0,nyu
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq(iy0,ix0,iz) &
                         +w10(iyf,ixf)*qq(iy1,ix0,iz) &
                         +w01(iyf,ixf)*qq(iy0,ix1,iz) &
                         +w11(iyf,ixf)*qq(iy1,ix1,iz)
    enddo
  enddo

   !Convert to real*4 (single precision):
  do ix=0,nxum1
    do iy=0,nyu
      qar4(iy,ix)=real(qa(iy,ix))
    enddo
  enddo

   !Write data for this layer:
  write(*,'(a,i1,a)') ' Writing output file for layer ',iz,' ...'
  write(44,rec=iz) real(t),qar4

enddo
 !Ends loop over layers

close(44)

write(*,*)
write(*,*) '   Image the file using '
write(*,'(a,i5,1x,i5,a)') ' dataview -ndim ',nxu,nyu+1,' fine/qq'//pind//'.r4'

return
end subroutine

!==============================================

end program
