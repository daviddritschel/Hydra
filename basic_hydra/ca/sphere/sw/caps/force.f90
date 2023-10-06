module force

!===========================================================================
! Module which generates a random field with a selected range of spherical
! harmonics (in degree n) with normalisation norm(n).

! Initialise with
  
! call init_forcing(norm)

! where norm is a double precision array dimensioned norm(nbeg:nend)

! Subsequently, generate a random field ff on a latitude-longitude grid by

! call generate_forcing(ff, frms)

! where frms is the r.m.s. value of ff, while ff is a double precision
! array dimensioned ff(ng,nt).

! David G. Dritschel, St Andrews UK, 16 December 2021
!===========================================================================

 !Import parameters and constants:
use constants

 !Import spherical module (stand alone):
use spherical

implicit none

 !Global variables:
double precision:: poly(ng,ntot), wave(ntot)
double precision:: sinphi(ng), cosphi(ng), lambda(nt)
double precision:: rmsfac

contains

!=========================================================================
  
subroutine init_forcing(norm)

implicit none

 !Passed variable:
double precision:: norm(nbeg:nend)

 !Local variables:
double precision:: p(ng), pfac

integer, dimension(:), allocatable :: seed
integer:: i, j, k, m, n

!-----------------------------------------------------------------
 !Define half-grid values of sinphi = sin(latitude):
do j = 1, ng
  sinphi(j) = sin(dl*(dble(j)-0.5d0)-hpi)
enddo

 !Define cos(lat):
cosphi = sqrt(1.d0-sinphi**2)

 !Needed to compute field r.m.s. value in generate_forcing:
rmsfac=1.d0/((f1112*(cosphi(1)+cosphi(ng))+sum(cosphi(2:ng-1)))*dble(nt))

 !Define longitudes:
do i = 1, nt
  lambda(i) = dl*dble(i-1)-pi
enddo

 !Initialise random number generator:
call random_seed(size = k)
allocate(seed(1:k))
seed(:) = iseed
do i = 1, iseed
  call random_seed(put = seed)
enddo

!-----------------------------------------------------------------
 !Generate normalised associated Legendre polynomials:
k = 0
do n = nbeg, nend
   !Normalisation factor:
  pfac = norm(n)

   !m = 0; include a factor of 1/2 (usual cosine series weight):
  k = k+1                         ! Index to store the polynomials
  wave(k) = 0.d0                  ! Stores azimuthal wavenumber
  call generate_spherical(sinphi,ng,n,0,p)
   !p = normalised Y_n^0 without the azimuthal dependence
  poly(:,k) = 0.5d0*pfac*p        ! Stores re-scaled polynomial

   !m > 0:
  do m = 1, n
    k = k+1
    wave(k) = dble(m)
    call generate_spherical(sinphi,ng,n,m,p)
     !p = normalised Y_n^m without the azimuthal dependence
    poly(:,k) = pfac*p
  enddo

enddo

return
end subroutine init_forcing

!=========================================================================

subroutine generate_forcing(ff, frms)
 !Generates a field ff with a narrow-band spectum and on a sphere.
 !The r.m.s. value of ff is passed in frms 

implicit none

 !Passed variables:
double precision:: ff(ng,nt), frms

 !Local variables:
double precision:: wka(ng,nt), cosfac(nt)
double precision:: theta, sfac
integer:: i, k

!-----------------------------------------------------------------
 !Generate random field:
ff = 0.d0
do k = 1,ntot
   !Random phase (theta):
  call random_number(theta)       !   0 < theta < 1 here
  theta = pi*(2.d0*theta-1.d0)    ! -pi < theta < pi now

   !cos(m*lambda + theta):
  cosfac = cos(wave(k)*lambda+theta)

   !Sum poly*cosfac into ff:
  do i = 1,nt
     ff(:,i) = ff(:,i)+poly(:,k)*cosfac(i)
  enddo
enddo

 !Scale ff to have correct rms value:
sfac = frms/rmsval(ff)
ff = sfac*ff

return
end subroutine generate_forcing

!=========================================================================

function rmsval(var)
! Computes the rms value of a gridded field var.

implicit none

 !Function name returns result:
double precision:: rmsval

 !Passed variable:
double precision:: var(ng,nt)

 !Local variables:
double precision:: wka(ng,nt)
integer:: i

 !-------------------------------------------------
do i = 1,nt
  wka(:,i) = cosphi*var(:,i)**2
enddo
 !cosphi is cos(latitude) here: the area weight.

rmsval = 0.d0
do i = 1,nt
  rmsval = rmsval+f1112*(wka(1,i)+wka(ng,i))+sum(wka(2:ng-1,i))
enddo
rmsval = sqrt(rmsfac*rmsval)
 !See subroutine init_forcing above for rmsval

end function rmsval

end module force
