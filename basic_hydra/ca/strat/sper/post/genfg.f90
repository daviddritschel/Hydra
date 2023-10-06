program genfg

! Program for post-processing the contour and residual files created by a run with 
! the casl suite of f90 code.

use contours
use congen

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Declarations:
real:: tr4,zdr4(0:ny,0:nxm1)
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

  do ix=0,nxm1
    do iy=0,ny
      bb(iy,ix)=zero
    enddo
  enddo

  call ugridsave(bb,xb,yb,bjump,bavg,nextb,indb,npb,i1b,nb,nptb,t,0)

else
   !Open residual needed to build ultra-fine-grid vorticity with congen:
  open(40,file='cont/zzresi.r4',form='unformatted', &
      & access='direct',status='old',recl=nbytes)
  read(40,rec=iind) tr4,zdr4
  close(40) 
  do ix=0,nxm1
    do iy=0,ny
      zd(iy,ix)=dble(zdr4(iy,ix))
    enddo
  enddo

   !Open buoyancy contours:
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

  call ugridsave(zd,xz,yz,zjump,zavg,nextz,indz,npz,i1z,nz,nptz,t,1)

endif

contains 

!-----------------------------------------------------------------
subroutine ugridsave(qq,xq,yq,dq,qavg,nextq,indq,npq,i1q,nq,nptq,t,ires)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
integer,parameter:: dbleint=selected_int_kind(16)
integer(kind=dbleint),parameter:: npu=int(nxu,kind=dbleint)*int(nyu+1,kind=dbleint)
real:: qar4(0:nyu,0:nxum1)
double precision:: qq(0:ny,0:nxm1)
double precision:: xq(npm),yq(npm)
integer:: nextq(npm),indq(nm),npq(nm),i1q(nm)

if (nq .gt. 0) call con2ugrid(xq,yq,dq,qavg,nextq,nptq)

if (ires .eq. 1) then
  !Bi-linear interpolate the residual qq to the fine grid and add to qa:
  do ix=0,nxum1
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

do ix=0,nxum1
  do iy=0,nyu
    qar4(iy,ix)=real(qa(iy,ix))
  enddo
enddo

write(*,*) ' Writing output file...'

if (ires .eq. 1) then
  open(44,file='fine/zz'//pind//'.r4',form='unformatted', &
      & access='direct',status='replace',recl=4+npu*4)
else
  open(44,file='fine/bb'//pind//'.r4',form='unformatted', &
      & access='direct',status='replace',recl=4+npu*4)
endif

write(44,rec=1) real(t),qar4
close(44)

write(*,'(a,i8,a,i8)') ' Fine grid dimensions are: ',nxu,', by ',nyu

return
end subroutine

!==============================================

end program
