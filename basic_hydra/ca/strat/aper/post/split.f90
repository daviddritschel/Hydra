program split

! Program for computing the division of kinetic energy and enstrophy
! into a large scale part and a small scale part. (uses subroutine filter)

use evolution
use spectral

implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: nbread=(nx+1)*(ny+1)+1
real:: tin,zzin(0:ny,0:nx)
double precision:: zzlo(0:ny,0:nx),zzhi(0:ny,0:nx)

!---------------------------------------------------------

 !Select frame to process:
write(*,*) ' Number of 1-2-1 filter passes to use:'
read(*,*) nrep

call init_spectral(one)

 !Open vorticity data file:
open(15,file='zz.r4',form='unformatted', & 
        & access='direct',status='old',recl=nbread*4)

!------code to find zzlo, zzhi for one frame------------
!     This code is older and would need adpation if it is 
!       to be used...!!!!
! !Select frame to process:
!write(*,*) ' Time frame to process (1,2,...)?:'
!read(*,*) kfr
!nrec=1+(nx+1)*(ny+1)
!
!if (kfr .gt. 1) then
!  do j=1,nrec*(kfr-1)
!    read(15,*)
!  enddo
!endif
!
!read(15,*) t
!do ix=0,nx
!  do iy=0,ny
!    read(15,*) zzlo(iy,ix)
!    zzhi(iy,ix)=zzlo(iy,ix)
!  enddo
!enddo
!close(15)
!
! !Low-pass filter zz:
!call filter(zzlo,0,nrep)
!
! !Calculate hi-pass filtered zz:
!do ix=0,nx
!  do iy=0,ny
!    zzhi(iy,ix)=zzhi(iy,ix)-zzlo(iy,ix)
!  enddo
!enddo
!
! !Write out filtered fields
!open(16,file='zzlo.dat',status='unknown')
!open(17,file='zzhi.dat',status='unknown')
!write(16,'(f12.5)') t
!write(17,'(f12.5)') t
!do ix=0,nx
!  do iy=0,ny
!    write(16,'(f16.7)') zzlo(iy,ix)
!    write(17,'(f16.7)') zzhi(iy,ix)  
!  enddo
!enddo
!
!close(16)
!close(17)
!-----------------------------------------------------------

open(17,file='split.asc',status='unknown')
loop=0
iread=0
 !Read field and process:
do 
  loop=loop+1
  iread=0
  read(15,rec=loop,iostat=iread) tin,zzin
  if (iread .ne. 0) exit 
  t=dble(tin)
  do ix=0,nx
    do iy=0,ny
      zzlo(iy,ix)=dble(zzin(iy,ix))
      zzhi(iy,ix)=zzlo(iy,ix)
    enddo
  enddo

   !Compute kinetic energy & enstrophy of full zz:
  call main_invert(zzlo,uu,vv)
  call l2norm(zzlo,zzl2)
  call l2norm(uu,uul2)
  call l2norm(vv,vvl2)
  eke=f12*(uul2+vvl2)

   !Low-pass filter zz:
  call filter(zzlo,0,nrep)

   !Calculate hi-pass filtered zz:
  do ix=0,nx
    do iy=0,ny
      zzhi(iy,ix)=zzhi(iy,ix)-zzlo(iy,ix)
    enddo
  enddo
   !Compute enstrophy split:
  call l2norm(zzlo,zzlol2)
  call l2norm(zzhi,zzhil2)

   !Compute kinetic energy of low pass filtered zz:
  call main_invert(zzlo,uu,vv)
  call l2norm(uu,uul2)
  call l2norm(vv,vvl2)
  ekelo=f12*(uul2+vvl2)

   !Compute kinetic energy of hi pass filtered zz:
  call main_invert(zzhi,uu,vv)
  call l2norm(uu,uul2)
  call l2norm(vv,vvl2)
  ekehi=f12*(uul2+vvl2)

   !Output energy split to a file and screen:
  eke=eke+small
  zzl2=zzl2+small
  write(17,'(f7.2,4(2x,f10.8))') t,ekelo/eke,ekehi/eke,zzlol2/zzl2,zzhil2/zzl2
  write( *,'(f7.2,4(2x,f10.8))') t,ekelo/eke,ekehi/eke,zzlol2/zzl2,zzhil2/zzl2
enddo
close(17)

end program
