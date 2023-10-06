program jets

! Counts the number of jets with u > u_rms

use parameters

implicit none

double precision,parameter:: frms=1.d0, frac=0.1d0
double precision,parameter:: f1112=11.d0/12.d0
double precision:: zuu(ng), phi(ng), clat(ng)
double precision:: tmp(ng), pj(ng)
double precision:: tc, t, oma, um, adpj, dpj, dpjm, anj, rsum, rsumi
integer:: iread, j, nj, jm, na

!---------------------------------------------------------------------
! Open all files:
open(51,file='evolution/zu.asc',status='old')
open(61,file='jets.asc',status='replace')

write(*,*) 'Time from which to compute time-average number of jets?'
read(*,*) tc
write(*,*)

! Initialise:
anj=0.d0
na=0
read(51,*) t
do j=1,ng
  read(51,*) zuu(j),phi(j)
enddo
rewind(51)
! Needed only once to define um below:
clat=cos(phi)
rsum=f1112*(clat(1)+clat(ng))+sum(clat(2:ng-1))
rsumi=1.d0/rsum
  
! Read data and process:
do
  iread=0
  read(51,*,iostat=iread) t
  if (iread .ne. 0) exit 
  ! Read zonal-average zonal velocity, zuu:
  do j=1,ng
    read(51,*) zuu(j)
  enddo

  ! Define minimum u_bar to be a jet:
  oma=(f1112*(zuu(1)+zuu(ng))+sum(zuu(2:ng-1)))*rsumi
  ! oma above is the mean angular velocity
  zuu=zuu-oma*clat
  ! We subtract oma*cos(phi) above to define residual zuu
  tmp=clat*zuu**2
  ! Compute rms value of zuu and multiply be frms (see parameters above):
  um=frms*sqrt((f1112*(tmp(1)+tmp(ng))+sum(tmp(2:ng-1)))*rsumi)+1.d-8

  ! Count jets:
  nj=0
  do j=2,ng-1
    if ((zuu(j)-zuu(j-1) > 0) .and. (zuu(j)-zuu(j+1) > 0)) then
      ! zuu(j) is a local maximum
      if (zuu(j) > um) then
        ! This is big enough to be a jet:
        nj=nj+1
        pj(nj)=phi(j)
      endif
    endif
  enddo

  ! Possibly merge close jets:
  if (nj > 1) then
    do
      dpjm=pi
      ! Find closest pair of jets:
      do j=1,nj-1
        dpj=pj(j+1)-pj(j)
        if (dpj < dpjm) then
          jm=j
          dpjm=dpj
        endif
      enddo

      ! Exit if this closest pair is larger than a certain fraction
      ! of the average jet spacing (adpj):
      adpj=pi/dble(nj+1)
      if (dpjm > frac*adpj) exit
      write(*,'(a,f8.1)') ' Merging a pair of jets at t = ',t
      ! Otherwise merge jets jm and jm+1:
      pj(jm)=0.5d0*(pj(jm)+pj(jm+1))
      ! Pack rest of the array so it is contiguous:
      if (jm < nj-1) then
        do j=jm+1,nj-1
          pj(j)=pj(j+1)
        enddo
      endif
      ! Reduce jet count:
      nj=nj-1
      if (nj == 1) exit
      ! Otherwise, go back and process next closest pair:
    enddo
  endif

  ! Write data:
  write(61,*) t,nj

  if (t+1.d-6 >= tc) then
    ! Accumulate time average:
    anj=anj+dble(nj)
    na=na+1
  endif

enddo

! Finalise average:
anj=anj/dble(na)

! Close files:
close(51)
close(61)

write(*,*)
write(*,'(a,f5.2)') ' The time-average number of jets is ',anj
write(*,*)
write(*,*) 'The number of jets vs t is in jets.asc'

! End main program
end program jets
