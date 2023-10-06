program fgrid

! Program for post-processing the contour and residual files created by a run with 
! the casl suite of f90 code.

use contours
use congen

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Declarations:
character*3:: pind
integer :: iclosed(nm)

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
  open(40,file='contours/bbcont'//pind,status='old')
  read(40,*) nb,nptb,t,bjump,bavg
  do j=1,nb
    read(40,*) npb(j),i1b(j),indb(j),iclosed(j)
   !iclosed = 0 for an open contour, and 1 for a closed one
  enddo
  do i=1,nptb
    read(40,*) xb(i),yb(i)
  enddo
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

  do ix=0,nx
    do iy=0,ny
      bb(iy,ix)=zero
    enddo
  enddo

  call ugridsave(bb,xb,yb,bjump,bavg,nextb,indb,npb,i1b,nb,nptb,t,0)

else
   !Open residual needed to build ultra-fine-grid vorticity with congen:
  open(40,file='contours/zzresi'//pind,status='unknown')
  read(40,*) t
  do ix=0,nx
    do iy=0,ny
      read(40,*) zd(iy,ix)
    enddo
  enddo
  close(40) 

   !Open vorticity contours:
  open(40,file='contours/zzcont'//pind,status='unknown')
  read(40,*) nz,nptz,t,zjump,zavg
  do j=1,nz
    read(40,*) npz(j),i1z(j),indz(j),iclosed(j)
   !iclosed = 0 for an open contour, and 1 for a closed one
  enddo
  do i=1,nptz
    read(40,*) xz(i),yz(i)
  enddo
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

  call ugridsave(zd,xz,yz,zjump,zavg,nextz,indz,npz,i1z,nz,nptz,t,1)
 
endif


contains 

!-----------------------------------------------------------------
subroutine ugridsave(qq,xq,yq,dq,qavg,nextq,indq,npq,i1q,nq,nptq,t,ires)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(0:ny,0:nx)
double precision:: xq(npm),yq(npm)
integer:: nextq(npm),indq(nm),npq(nm),i1q(nm)

call con2ugrid(xq,yq,dq,qavg,nextq,nptq)

if (ires .eq. 1) then
  !Bi-linear interpolate the residual qq to the fine grid and add to qa:
  do ix=0,nxu
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    do iy=0,nyu
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq(iy0,ix0) &
                       & +w10(iyf,ixf)*qq(iy1,ix0) &
                       & +w01(iyf,ixf)*qq(iy0,ix1) &
                       & +w11(iyf,ixf)*qq(iy1,ix1)
    enddo
  enddo
endif

write(*,*) ' Writing output file...'
open(44,file='fine-grid.dat',status='unknown')
write(44,'(f8.2)') tgrid
do ix=0,nxu
  do iy=0,nyu
    write(44,'(f16.7)') qa(iy,ix)
  enddo
enddo
close(44)

return
end subroutine

!==============================================

end program
