program rest

!=============================================================
!  Sets up q = planetary vorticity, with topographic forcing
!=============================================================

 !Import spectral module:
use spectral

 !Declarations:
implicit none

double precision:: qq(ng,nt)
integer:: i

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

 !Define PV for a rest state:
do i=1,nt
  qq(:,i)=fpole*slat
enddo

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

 ! If topographic forcing is present, generate and save initial topography:
if (forcing) then
  call generate_forcing(qq,brms)
  open(11,file='bb_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,qq
  close(11)
endif

end program rest
