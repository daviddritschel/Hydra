program decomp
!  -------------------------------------------------------------------------
!  |   Decomposes the energy and enstrophy into the parts associated with  |
!  |   each (vertical) mode.                                               |
!  |                                                                       |
!  |   Output is to the formatted file "diagnostics/ez.asc" listing data   |
!  |   in the form                                                         |
!  |             t, K_1, K_2, K_tot, APE, E, Z_1, Z_2, Z                   |
!  |   with each record corresponding to a given time read from the input  |
!  |   files qq1.r4 & qq2.r4, and with K = kinetic energy, E = kinetic +   |
!  |   potential energy, APE = available potential energy, and             |
!  |   Z = enstrophy.                                                      |
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
double precision:: u1(0:ny,0:nxm1), v1(0:ny,0:nxm1), q1(0:ny,0:nxm1)
double precision:: u2(0:ny,0:nxm1), v2(0:ny,0:nxm1), q2(0:ny,0:nxm1)

 !Diagnostic quantities:
double precision:: k1, k2, ktot, ape, etot, z1, z2, ztot
double precision:: t, ul2, vl2, pl1, pl2, ql2
double precision:: qq1, qq2
real:: q1r4(0:ny,0:nxm1), q2r4(0:ny,0:nxm1), tr4

integer:: loop, iread, ix, iy

!---------------------------------------------------------------
 !Open files containing PV field in each layer:
open(31,file='evolution/qq1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(32,file='evolution/qq2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output file:
open(22,file='diagnostics/ez.asc',status='replace')

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

   !Compute mode 1 and mode 2 parts of u, v & q - beta*y:
  do ix=0,nxm1
    do iy=0,ny
      u1(iy,ix)=vec11*uu(iy,ix,1)+vec12*uu(iy,ix,2)
      u2(iy,ix)=vec21*uu(iy,ix,1)+vec22*uu(iy,ix,2)

      v1(iy,ix)=vec11*vv(iy,ix,1)+vec12*vv(iy,ix,2)
      v2(iy,ix)=vec21*vv(iy,ix,1)+vec22*vv(iy,ix,2)

       !Add topography (f_0*H_b/H_1) to PV in lower layer:
      qq1=qq(iy,ix,1)-bety(iy)+fhb(iy,ix)
      qq2=qq(iy,ix,2)-bety(iy)
      q1(iy,ix)=vec11*qq1+vec12*qq2
      q2(iy,ix)=vec21*qq1+vec22*qq2

       !Redefine lower-layer streamfunction for use below:
      pp(iy,ix,1)=pp(iy,ix,1)-alpha*pp(iy,ix,2)
    enddo
  enddo

   !Compute available potential energy from interface displacements:
  call l2norm(pp(0,0,1),pl1)
  call l2norm(pp(0,0,2),pl2)
  ape=f12*h1h2kdbarsq*(pl1+alpha*alphac*pl2)
   !alphac = 1 - alpha here

   !Compute kinetic energy and enstrophy in each mode:
  call l2norm(u1,ul2)
  call l2norm(v1,vl2)
  k1=coeff1*(ul2+vl2)

  call l2norm(u2,ul2)
  call l2norm(v2,vl2)
  k2=coeff2*(ul2+vl2)
  ktot=k1+k2
  etot=ktot+ape

  call l2norm(q1,ql2)
  z1=coeff1*ql2

  call l2norm(q2,ql2)
  z2=coeff2*ql2
  ztot=z1+z2

   !Write diagnostic data for this time:
  write(22,'(f8.2,8(1x,f11.7))') t, k1, k2, ktot, ape, etot, z1, z2, ztot
enddo

 !Close output file:
close(22)

write(*,*)
write(*,*) ' The results are ready in ez.asc; this contains a listing of'
write(*,*) '       t, K_1, K_2, K, APE, E, Z_1, Z_2, and Z'
write(*,*) ' where K = kinetic energy, E = kinetic + potential energy, '
write(*,*) ' APE = available potential energy, and Z = (potential) entrophy.'
write(*,*) ' NOTE: the subscripts 1 or 2 refer to the vertical mode.'
write(*,*)
write(*,*) ' Plot using the command:  plotcol ez.asc'

return
end subroutine

 !End main program
end program
!=======================================================================
