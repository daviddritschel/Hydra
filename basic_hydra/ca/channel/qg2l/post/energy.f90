program energy
!  -------------------------------------------------------------------------
!  |   Calculates various energy components (kinetic & potential, zonal &  |
!  |   non-zonal, mode 1 & mode 2).                                        |
!  |                                                                       |
!  |   Output is to the formatted files "kpzonal.asc" & "kpeddy.asc"       |
!  |   listing data for each time t in the form                            |
!  |             t, K_1, K_2, P                                            |
!  |   where K_m is the kinetic energy in mode m and P is the available    |
!  |   potential energy.  kpzonal lists the energies due to the zonal mean |
!  |   streamfunction and velocity, while kpeddy lists those due to the    |
!  |   departures of streamfunction and velocity from their zonal mean.    |
!  -------------------------------------------------------------------------

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use spectral
use generic

implicit none
double precision:: fhb(0:ny,0:nxm1)
double precision:: dum
integer:: ix,iy

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
  open(12,file='topo.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
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
double precision:: u1bar(0:ny),u2bar(0:ny)
double precision:: p1bar(0:ny),p2bar(0:ny)

 !Diagnostic quantities:
double precision:: t, zfac, aac, zke1, zke2, zpe, eke1, eke2, epe
real:: q1r4(0:ny,0:nxm1), q2r4(0:ny,0:nxm1), tr4

integer:: loop, iread, ix, iy

!---------------------------------------------------------------
 !Open files containing PV field in each layer:
open(31,file='evolution/qq1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(32,file='evolution/qq2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output files:
open(21,file='diagnostics/kpzonal.asc',status='replace')
open(22,file='diagnostics/kpeddy.asc' ,status='replace')

zfac=one/dble(nx)
aac=alpha*alphac
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

  do ix=0,nxm1
    do iy=0,ny
      qq(iy,ix,1)=dble(q1r4(iy,ix))
      qq(iy,ix,2)=dble(q2r4(iy,ix))
    enddo
  enddo

   !Invert PV to get velocity field (uu,vv) and streamfunction pp:
  call main_invert(qq,fhb,t,uu,vv,pp)

   !Initialise zonal mean values:
  do iy=0,ny
    u1bar(iy)=zero
    u2bar(iy)=zero
    p1bar(iy)=zero
    p2bar(iy)=zero
  enddo

   !Compute mode 1 and mode 2 parts of u & v:
  do ix=0,nxm1
    do iy=0,ny
      u1(iy,ix)=vec11*uu(iy,ix,1)+vec12*uu(iy,ix,2)
      u2(iy,ix)=vec21*uu(iy,ix,1)+vec22*uu(iy,ix,2)
      u1bar(iy)=u1bar(iy)+u1(iy,ix)
      u2bar(iy)=u2bar(iy)+u2(iy,ix)

      v1(iy,ix)=vec11*vv(iy,ix,1)+vec12*vv(iy,ix,2)
      v2(iy,ix)=vec21*vv(iy,ix,1)+vec22*vv(iy,ix,2)

       !Redefine lower-layer streamfunction for use below in potential energy:
      pp(iy,ix,1)=pp(iy,ix,1)-alpha*pp(iy,ix,2)
      p1bar(iy)=p1bar(iy)+pp(iy,ix,1)
      p2bar(iy)=p2bar(iy)+pp(iy,ix,2)
    enddo
  enddo

   !Finalise zonal mean values:
  do iy=0,ny
    u1bar(iy)=zfac*u1bar(iy)
    u2bar(iy)=zfac*u2bar(iy)
    p1bar(iy)=zfac*p1bar(iy)
    p2bar(iy)=zfac*p2bar(iy)
  enddo

   !Compute zonal energies:
  zke1=f12*(u1bar(0)**2+u1bar(ny)**2)
  zke2=f12*(u2bar(0)**2+u2bar(ny)**2)
  zpe =f12*(p1bar(0)**2+p1bar(ny)**2 + aac*(p2bar(0)**2+p2bar(ny)**2))
  do iy=1,nym1
    zke1=zke1+u1bar(iy)**2
    zke2=zke2+u2bar(iy)**2
    zpe =zpe +p1bar(iy)**2 + aac*p2bar(iy)**2
  enddo
  zke1=f12*coeff1*ellx*gly*zke1
  zke2=f12*coeff2*ellx*gly*zke2
  zpe =f12*h1h2kdbarsq*ellx*gly*zpe
   !coeff1 & coeff2 ensure that the sum of the modal energies equals the 
   !sum of the (mass and thickness weighted) layer energies.  The factor
   !ellx is needed so that energies are area integrals (not domain averages).

   !Compute zonal energies after removing zonal means:
  do ix=0,nxm1
    do iy=0,ny
      u1(iy,ix)=u1(iy,ix)-u1bar(iy)
      u2(iy,ix)=u2(iy,ix)-u2bar(iy)
      pp(iy,ix,1)=pp(iy,ix,1)-p1bar(iy)
      pp(iy,ix,2)=pp(iy,ix,2)-p2bar(iy)
    enddo
  enddo

  eke1=zero
  eke2=zero
  epe =zero
  do ix=0,nxm1
    eke1=eke1+u1(0,ix)**2+u1(ny,ix)**2 + v1(0,ix)**2+v1(ny,ix)**2
    eke2=eke2+u2(0,ix)**2+u2(ny,ix)**2 + v2(0,ix)**2+v2(ny,ix)**2
    epe =epe +pp(0,ix,1)**2+pp(ny,ix,1)**2 + aac*(pp(0,ix,2)**2+pp(ny,ix,2)**2)
  enddo
  eke1=f12*eke1
  eke2=f12*eke2
  epe =f12*epe

  do ix=0,nxm1
    do iy=1,nym1
      eke1=eke1+u1(iy,ix)**2 + v1(iy,ix)**2
      eke2=eke2+u2(iy,ix)**2 + v2(iy,ix)**2
      epe =epe +pp(iy,ix,1)**2 + aac*pp(iy,ix,2)**2
    enddo
  enddo

  eke1=f12*coeff1*garea*eke1
  eke2=f12*coeff2*garea*eke2
  epe =f12*h1h2kdbarsq*garea*epe
   !Note: garea is the grid box area, glx*gly

   !Write diagnostic data for this time:
  write(21,'(f8.2,3(1x,f13.8))') t, zke1, zke2, zpe
  write(22,'(f8.2,3(1x,f13.8))') t, eke1, eke2, epe
enddo

 !Close files:
close(21)
close(22)
close(31)
close(32)

write(*,*)
write(*,*) '    The zonal modal kinetic energies and APE are in kpzonal.asc'
write(*,*) ' and the eddy modal kinetic energies and APE are in kpeddy.asc'

return
end subroutine

 !End main program
end program
!=======================================================================
