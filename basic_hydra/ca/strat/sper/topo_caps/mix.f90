program mix

! Computes the mean physical y coordinate (Y_bar) as a function of 
! buoyancy b from data in bb.r4, previously generated from topo_caps.
! Also computes the standard deviation sqrt<(Y-Y_bar)^2> as well as
! the bouyancy pdf, P(b).

use spectral

implicit none

 !Declarations:
double precision:: hh(0:nxm1),range
real bb(0:ny,0:nxm1)
real pdf(ncontb),bbar(ncontb),ybar(ncontb),yrms(ncontb)
real bbmin,bbmax,bjump,dbi,t,tbeg,tend,pdfsum,fnorm
integer:: ib,iread,ix,iy,k,kbeg,kend

!-----------------------------------------------------------------
 !Open file containing the buoyancy:
open(20,file='bb.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

 !Compute contour interval (bjump) for buoyancy from initial data:
read(20,rec=1,iostat=iread) t,bb
bbmax=0.
bbmin=0.
do ix=0,nxm1
  do iy=0,ny
    bbmax=max(bbmax,bb(iy,ix))
    bbmin=min(bbmin,bb(iy,ix))
  enddo
enddo
write(*,'(a,f9.5)') ' Minimum b used for diagnostics = ',bbmin
write(*,'(a,f9.5)') ' Maximum b  "    "       "      = ',bbmax
write(*,*)

 !Read bottom topography:
open(13,file='hh.asc',status='old')
do ix=0,nxm1
  read(13,*) hh(ix)
enddo
close(13)

 !Initialise conformal factor and coordinates:
range=bbmax-bbmin
call init_spectral(range,hh)

 !Adjust bbmin & bbmax slightly to fit data fully within the bins:
bbmin=bbmin-small*range
bbmax=bbmax+small*range
bjump=(bbmax-bbmin)/dble(ncontb)
dbi=1./bjump

write(*,*) ' Time interval to use for averaging?'
read(*,*) tbeg,tend
kbeg=nint(tbeg/tgsave)+1
kend=nint(tend/tgsave)+1

! Initialise bins for forming pdf, etc:
do ib=1,ncontb
  pdf(ib)=0.
  ybar(ib)=0.
  yrms(ib)=0.
enddo

!---------------------------------------------------------------
 !Read data and process:
do k=kbeg,kend
  read(20,rec=k,iostat=iread) t,bb
  if (iread .eq. 0) then
    write(*,'(a,f13.5)') ' Processing t = ',t

    do ix=0,nxm1
      do iy=0,ny
        ib=max(1,min(ncontb,1+int(dbi*(bb(iy,ix)-bbmin))))
        pdf(ib)=pdf(ib)+confac(iy,ix)
        ybar(ib)=ybar(ib)+confac(iy,ix)*yori(iy,ix)
        yrms(ib)=yrms(ib)+confac(iy,ix)*yori(iy,ix)**2
      enddo
    enddo
     !Note: the area of a grid box is proportional to confac(iy,ix)

  endif
enddo
close(20)

 !Get Y_bar and Y_rms, and normalise pdf so that int{P*db} = 1:
pdfsum=0.
do ib=1,ncontb
  bbar(ib)=bbmin+bjump*(float(ib)-0.5)
  ybar(ib)=ybar(ib)/pdf(ib)
  yrms(ib)=sqrt(abs(yrms(ib)/pdf(ib)-ybar(ib)**2))
  pdfsum=pdfsum+pdf(ib)
enddo

fnorm=1./(pdfsum*bjump)
do ib=1,ncontb
  pdf(ib)=pdf(ib)*fnorm
enddo

 !Write data:
open(30,file='Y.asc',status='replace')
open(40,file='P.asc',status='replace')
do ib=1,ncontb
  if (pdf(ib) .gt. 0.) then
    write(30,'(3(1x,f11.7))') bbar(ib),ybar(ib),yrms(ib)
    write(40,'(1x,f11.7,1x,f13.7)') bbar(ib),pdf(ib)
  endif
enddo
close(30)
close(40)

write(*,*) ' b vs Y_bar(b) & Y_rms(b) are listed in Y.asc, while'
write(*,*) ' b vs P(b) is listed in P.asc'

end program

