program ranpv

!=====================================================
!    Sets up a random PV relative vorticity field
!=====================================================

 !Import spectral module:
use spectral

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: zz(ng,nt)

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

 !Initialize random # generator:
do i=1,iseed
  uni=rand(0)
enddo

 !--------------------------------------------------------
write(*,*) ' Enter the rms relative vorticity:'
read(*,*) eps

 !generate random relative vorticity field zz with this rms value:
call ranspec(zz,eps)

 !Return zz to physical space:
call revfft(ng,nt,zz,trig,factors)

if (omega .ne. zero) then
   !Add f:
  do i=1,nt
    do j=1,ng
      zz(j,i)=zz(j,i)+cof(j)
    enddo
  enddo
endif

 !Write initial absolute vorticity field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

end program
