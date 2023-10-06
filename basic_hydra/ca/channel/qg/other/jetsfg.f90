program monfg
!  --------------------------------------------------------------------------
!  |   Takes the fine grid PV field generated by genfg.f90 and creates a    |
!  |   monotonic (in y) field (for each x) by sorting.                      |
!  |   Writes fine/monqq<pind>.r4 where <pind> is the selected time period  |
!  --------------------------------------------------------------------------

 !Import constants:
use constants

implicit none

 !Declarations:
integer,parameter:: dbleint=selected_int_kind(16)
integer(kind=dbleint),parameter:: npu=int(nxu,kind=dbleint)*int(nyu+1,kind=dbleint)

real:: tr4, qar4(0:nyu,0:nxum1)
real:: qmin, qmax, dbin, dbini
real:: cmult, sumj, stair(nz)
real:: pdf(nyu)

integer:: iind, ix, iy, iz, k
integer:: njets(nz)

logical struct
character(len=3):: pind

!-----------------------------------------------------------------
write(*,*) ' What is the numbered suffix of the fine-grid PV file?'
read(*,*) iind
write(pind(1:3),'(i3.3)') iind

 !Open input file:
open(33,file='fine/qq'//pind//'.r4',form='unformatted', &
    & access='direct',status='old',recl=4+npu*4)

write(*,*) ' Threshold probability as a multiple C of the mean probability?'
read(*,*) cmult

 !Process each layer in turn:
do iz=1,nz
   !Read PV field for this layer:
  read(33,rec=iz) tr4,qar4

   !Compute min & max PV:
  qmin=qar4(0,0)
  qmax=qar4(0,0)
  do ix=0,nxum1
    do iy=1,nyu
      qmin=min(qmin,qar4(iy,ix))
      qmax=max(qmax,qar4(iy,ix))
    enddo
  enddo

   !Increment in PV used for binning PV in each column ix = constant:
  dbin=(qmax-qmin)/float(nyu-1)
  dbini=1./dbin

   !Initialise probabilities:
  do k=1,nyu
    pdf(k)=0.
  enddo

  do ix=0,nxum1
     !Edge values at iy = 0 & nyu count half:
    k=nint((qar4(0  ,ix)-qmin)*dbini)+1
    pdf(k)=pdf(k)+0.5
    k=nint((qar4(nyu,ix)-qmin)*dbini)+1
    pdf(k)=pdf(k)+0.5

     !Interior values count 1:
    do iy=1,nyu-1
      k=nint((qar4(iy,ix)-qmin)*dbini)+1
      pdf(k)=pdf(k)+1.
    enddo
  enddo

   !Normalise probabilities so they integrate to 1 over a rescaled
   !q range going from 0 to 1 (without loss of generality):
  sumj=1.d-12
  do k=1,nyu
    sumj=sumj+pdf(k)
  enddo
  sumj=float(nyu)/sumj
  do k=1,nyu
    pdf(k)=pdf(k)*sumj
  enddo

   !Note, the mean probability is 1 by construction.

   !Diagnose layer:
  njets(iz)=0
  stair(iz)=0.
  struct=.false.
  do k=1,nyu
    if (pdf(k) .gt. cmult) then
      stair(iz)=stair(iz)+pdf(k)-cmult
      if (.not. struct) then
        njets(iz)=njets(iz)+1
        struct=.true.
      endif
    else
      struct=.false.
    endif
  enddo
  if (njets(iz) .gt. 0) njets(iz)=njets(iz)-1
  stair(iz)=stair(iz)/float(nyu)

enddo

write(*,*)
write(*,'(a,i3)') ' Number of jets in layer 1 = ',njets(1)
write(*,'(a,f9.7)')     ' Degree of staircasing = ',stair(1)

write(*,*)
write(*,'(a,i3)') ' Number of jets in layer 2 = ',njets(2)
write(*,'(a,f9.7)')     ' Degree of staircasing = ',stair(2)

close(33)

 !End main program
end program
!=======================================================================