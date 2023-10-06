program genXYfg

! From the topographic coordinates in coords.r8, this creates a fine
! grid version in the fine subdirectory suitable for use with topoview.

use contours

implicit none

 !Declarations:
integer,parameter:: dbleint=selected_int_kind(16)
integer(kind=dbleint),parameter:: npu=int(nxu+1,kind=dbleint)*int(nyu+1,kind=dbleint)
double precision:: xori(0:ny,0:nxm1),yori(0:ny,0:nxm1)
double precision:: xoriu(0:nyu,0:nxum1),yoriu(0:nyu,0:nxum1)
double precision:: dum
integer:: ix,iy,ixf,ix0,ix1,iyf,iy0,iy1

!-----------------------------------------------------------------
! Initialise weights for interpolation to the ultra-fine grid:
call init_contours

! Read coordinates of inversion grid:
open(20,file='coords.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(20,rec=1) dum,xori
read(20,rec=2) dum,yori
close(20)

! Interpolate X-x and Y-y:
do ix=0,nxm1
  do iy=0,ny
    xori(iy,ix)=xori(iy,ix)-(glx*dble(ix)-hlx)
    yori(iy,ix)=yori(iy,ix)-(gly*dble(iy)+ymin)
  enddo
enddo

do ix=0,nxum1
  ixf=ixfw(ix)
  ix0=ix0w(ix)
  ix1=ix1w(ix)

  do iy=0,nyu
    iyf=iyfw(iy)
    iy0=iy0w(iy)
    iy1=iy1w(iy)

    xoriu(iy,ix)=w00(iyf,ixf)*xori(iy0,ix0)+w10(iyf,ixf)*xori(iy1,ix0) &
              & +w01(iyf,ixf)*xori(iy0,ix1)+w11(iyf,ixf)*xori(iy1,ix1)

    yoriu(iy,ix)=w00(iyf,ixf)*yori(iy0,ix0)+w10(iyf,ixf)*yori(iy1,ix0) &
              & +w01(iyf,ixf)*yori(iy0,ix1)+w11(iyf,ixf)*yori(iy1,ix1)
 
  enddo
enddo

! Add back x & y:
do ix=0,nxum1
  do iy=0,nyu
    xoriu(iy,ix)=xoriu(iy,ix)+xgu(ix)
    yoriu(iy,ix)=yoriu(iy,ix)+ygu(iy)
  enddo
enddo

! Write coordinates of ultra-fine grid:
open(20,file='fine/coords.r8',form='unformatted', &
    & access='direct',status='replace',recl=8+npu*8)
write(20,rec=1) 0.d0,xoriu
write(20,rec=2) 0.d0,yoriu
close(20)

end program
