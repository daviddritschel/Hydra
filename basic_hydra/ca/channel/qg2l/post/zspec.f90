program zspec
!  -------------------------------------------------------------------------
!  |   Diagnoses the y-integrated zonal kinetic energy spectrum in each    |
!  |   vertical mode (1 or 2) as well as the y-integrated potential energy |
!  |   spectrum.                                                           |
!  |                                                                       |
!  |   Output is to formatted file "zspec.asc" which may be plotted using  |
!  |   zsv in the scripts subdirectory.                                    |
!  -------------------------------------------------------------------------

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use spectral

implicit none
double precision:: fhb(0:ny,0:nxm1)
double precision:: dum
integer:: ix,iy

!-----------------------------------------------------------------
 !Initialise inversion constants and arrays (module spectral.f90):
call init_spectral

 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
  open(12,file='topo.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) dum,fhb
  close(12)
else
  do ix=0,nxm1
    do iy=0,ny
      fhb(iy,ix)=zero
    enddo
  enddo
endif

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(0:ny,0:nxm1,nz),pp(0:ny,0:nxm1,nz)
double precision:: uu(0:ny,0:nxm1,nz),vv(0:ny,0:nxm1,nz)
double precision:: u1(0:ny,0:nxm1),v1(0:ny,0:nxm1)
double precision:: u2(0:ny,0:nxm1),v2(0:ny,0:nxm1)
double precision:: p1(0:ny,0:nxm1),p2(0:ny,0:nxm1)

 !Spectra:
double precision:: k1spec(0:nwx),k1spectavg(0:nwx)
double precision:: k2spec(0:nwx),k2spectavg(0:nwx)
double precision:: pespec(0:nwx),pespectavg(0:nwx)
double precision:: sk1spec,sk2spec,spespec

 !Single precision input data:
real:: q1r4(0:ny,0:nxm1),q2r4(0:ny,0:nxm1)
real:: tr4,t1,t2,avgfac

 !Local variables:
double precision:: trig(2*nx),rk(nx),alkx(nwx),scx,t
double precision:: ymult(0:ny)
double precision, parameter:: zfac=one/dble(ny)
integer:: factors(5)
integer:: loop,iread,ix,iy,kx,kc,ntavg
logical:: timeavg

!---------------------------------------------------------------
 !Open files containing PV field in each layer:
open(31,file='evolution/qq1.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
open(32,file='evolution/qq2.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

write(*,*) ' Time period to carry out time average [t1,t2]?'
read(*,*) t1,t2

timeavg=t2 .ge. t1
if (timeavg) then
   !Initialise time-averaged spectra:
  do kx=0,nwx
    k1spectavg(kx)=zero
    k2spectavg(kx)=zero
    pespectavg(kx)=zero
  enddo
  ntavg=0
  t1=t1-0.001d0*tgsave
  t2=t2+0.001d0*tgsave
endif

 !Open output file (view as an animation using zsv in scripts):
open(51,file='spectra/zspec.asc',status='replace')

 !Initialise 1D FFTs:
call initfft(nx,factors,trig)

 !Define log of x wavenumbers for output:
scx=twopi/ellx
do kx=1,nwx
  alkx(kx)=log10(scx*dble(kx))
enddo

 !Spectral averaging factors:
ymult(0)=f12*zfac
do iy=1,nym1
  ymult(iy)=zfac
enddo
ymult(ny)=f12*zfac

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,q1r4
  if (iread .ne. 0) exit 
  read(32,rec=loop,iostat=iread) tr4,q2r4
  t=dble(tr4)
  write(*,'(a,f13.5)') ' Processing t = ',t

   !Convert data to double precision:
  do ix=0,nxm1
    do iy=0,ny
      qq(iy,ix,1)=dble(q1r4(iy,ix))
      qq(iy,ix,2)=dble(q2r4(iy,ix))
    enddo
  enddo

   !Invert PV to get velocity field (uu,vv) and streamfunction (pp):
  call main_invert(qq,fhb,t,uu,vv,pp)

   !Compute mode 1 and mode 2 velocity components:
  do ix=0,nxm1
    do iy=0,ny
      u1(iy,ix)=vec11*uu(iy,ix,1)+vec12*uu(iy,ix,2)
      u2(iy,ix)=vec21*uu(iy,ix,1)+vec22*uu(iy,ix,2)
      v1(iy,ix)=vec11*vv(iy,ix,1)+vec12*vv(iy,ix,2)
      v2(iy,ix)=vec21*vv(iy,ix,1)+vec22*vv(iy,ix,2)
       !Used for obtaining APE from interface displacements:
      p1(iy,ix)=pp(iy,ix,1)-alpha*pp(iy,ix,2)
      p2(iy,ix)=pp(iy,ix,2)
    enddo
  enddo

   !FFT in x only:
  call forfft(nyp1,nx,u1(0,0),trig,factors) 
  call forfft(nyp1,nx,u2(0,0),trig,factors) 
  call forfft(nyp1,nx,v1(0,0),trig,factors) 
  call forfft(nyp1,nx,v2(0,0),trig,factors) 
  call forfft(nyp1,nx,p1(0,0),trig,factors) 
  call forfft(nyp1,nx,p2(0,0),trig,factors) 

   !Initialise spectra:
  do kx=0,nwx
    k1spec(kx)=zero
    k2spec(kx)=zero
    pespec(kx)=zero
  enddo

   !Average spectra over y:
  kx=0 !Zero wavenumber
  do iy=0,ny
    k1spec(kx)=k1spec(kx)+ymult(iy)*(u1(iy,kx)**2+v1(iy,kx)**2)
    k2spec(kx)=k2spec(kx)+ymult(iy)*(u2(iy,kx)**2+v2(iy,kx)**2)
    pespec(kx)=pespec(kx)+ymult(iy)*(p1(iy,kx)**2+alpha*alphac*p2(iy,kx)**2)
  enddo

  kx=nwx !Nyquist wavenumber
  do iy=0,ny
    k1spec(kx)=k1spec(kx)+ymult(iy)*(u1(iy,kx)**2+v1(iy,kx)**2)
    k2spec(kx)=k2spec(kx)+ymult(iy)*(u2(iy,kx)**2+v2(iy,kx)**2)
    pespec(kx)=pespec(kx)+ymult(iy)*(p1(iy,kx)**2+alpha*alphac*p2(iy,kx)**2)
  enddo

  do kx=1,nwx-1 !Add contributions from both cosines (kx) & sines (kc)
    kc=nx-kx
    do iy=0,ny
      k1spec(kx)=k1spec(kx)+ymult(iy)*(u1(iy,kx)**2+v1(iy,kx)**2 &
                                  &   +u1(iy,kc)**2+v1(iy,kc)**2)
      k2spec(kx)=k2spec(kx)+ymult(iy)*(u2(iy,kx)**2+v2(iy,kx)**2 &
                                  &   +u2(iy,kc)**2+v2(iy,kc)**2)
      pespec(kx)=pespec(kx)+ymult(iy)*( &
                     &    p1(iy,kx)**2+alpha*alphac*p2(iy,kx)**2 &
                     &   +p1(iy,kc)**2+alpha*alphac*p2(iy,kc)**2)
    enddo
  enddo

   !Re-scale the spectra and obtain sums for information only:
  sk1spec=zero
  sk2spec=zero
  spespec=zero
  do kx=0,nwx
    k1spec(kx)=coeff1*k1spec(kx)
    k2spec(kx)=coeff2*k2spec(kx)
    pespec(kx)=h1h2kdbarsq*pespec(kx)
    sk1spec=sk1spec+k1spec(kx)
    sk2spec=sk2spec+k2spec(kx)
    spespec=spespec+pespec(kx)
  enddo

  if (timeavg) then
     !Accumulate time average if t1-eps < t < t2+eps where eps = tgsave/1000:
    if (t1 .lt. t .and. t .lt. t2) then
      ntavg=ntavg+1
      do kx=0,nwx
        k1spectavg(kx)=k1spectavg(kx)+k1spec(kx)
        k2spectavg(kx)=k2spectavg(kx)+k2spec(kx)
        pespectavg(kx)=pespectavg(kx)+pespec(kx)
      enddo
    endif
  endif

   !Write out spectrum (only for kx > 0) to file:
  write(51,'(f8.2,3(1x,e14.7),1x,i5)') tr4,sk1spec,sk2spec,spespec,nwx
  do kx=1,nwx
    write(51,'(4(1x,f12.8))') alkx(kx),log10(k1spec(kx)+1.d-50), &
                               &       log10(k2spec(kx)+1.d-50), &
                               &       log10(pespec(kx)+1.d-50)
  enddo
   !Note: alk(kx) = log_10(kx)

enddo

if (timeavg) then
   !Finalise time-averaged spectra:
  avgfac=one/dble(ntavg)
  do kx=0,nwx
    k1spectavg(kx)=k1spectavg(kx)*avgfac
    k2spectavg(kx)=k2spectavg(kx)*avgfac
    pespectavg(kx)=pespectavg(kx)*avgfac
  enddo

   !Write data for use by plotcol, with y in the second row:
  open(41,file='spectra/k1zspectavg.asc',status='replace')
  open(42,file='spectra/k2zspectavg.asc',status='replace')
  open(43,file='spectra/pezspectavg.asc',status='replace')
  do kx=1,nwx
    write(41,'(2(1x,f12.8))') alkx(kx),log10(k1spectavg(kx))
    write(42,'(2(1x,f12.8))') alkx(kx),log10(k2spectavg(kx))
    write(43,'(2(1x,f12.8))') alkx(kx),log10(pespectavg(kx))
  enddo
   !Note: alk(kx) = log_10(kx)
  close(41)
  close(42)
  close(43)
endif

 !Close output file:
close(51)

write(*,*)
write(*,*) ' The results are ready in zspec.asc; for imaging, use'
write(*,*) ' zsv zspec.asc'
write(*,*) ' where zsv is in the scripts subdirectory.'

if (timeavg) then
  write(*,*)
  write(*,*) ' Time averaged KE_1, KE_2 & PE spectra are in k1zspectavg.asc,'
  write(*,*) ' k2zspectavg.asc & pezspectavg.asc.  Plot using'
  write(*,*)
  write(*,*) ' plotcol k1zspectavg.asc k2zspectavg.asc pezspectavg.asc'
endif

return
end subroutine

 !End main program
end program
!=======================================================================
