module spherical

!===========================================================================
! Module used to generate the normalised associated Legendre polynomials
! lambda_m^n appearing in the definition of spherical harmonics, i.e.
! Y_n^m = lambda_m^n * exp(i*m*theta) where theta is longitude.

! Generate lambda_n^m for an array of z grid points of length n using
!
! call generate_spherical(z, nz, n, m, p)
!
! where p = lambda_n^m evaluated at all z grid points.

! *** Both z and p are double precision arrays of length nz ***

! Developed from the Julia code by Justin Willmert found at
!https://justinwillmert.com/articles/2020/pre-normalizing-legendre-polynomials/
  
! David G. Dritschel, St Andrews UK, 13 December 2021
!===========================================================================

implicit none

contains

subroutine generate_spherical(z, nz, n, m, p)

implicit none

 !Passed variables:
integer:: nz, n, m
double precision:: z(nz), p(nz)

 !Local constants:
double precision,parameter:: pi = 3.141592653589793238462643383279502884197169399375105820974944592307816d0
double precision,parameter:: p00 = 1.d0/sqrt(4.d0*pi) ! Spherical normalisation

 !Local variables:
double precision:: r(nz),pmin1(nz),pmin2(nz)
double precision:: dn,fac
double precision:: mu,nu,alp,bet
integer:: mp,np,msq,nm1

!----------------------------------------------------------------------
! lambda_0^0:
p = p00
if (n == 0) return

! Increase index to m:
r = sqrt(1.d0-z**2)
do mp = 1, m
  mu = sqrt(1.d0+1.d0/dble(2*mp))
  p = -mu*r*p
enddo
if (n == m) return

nu = sqrt(dble(1+2*(m+1)))
pmin2 = p
pmin1 = nu*z*pmin2
if (n == m+1) then
  p = pmin1
  return
endif

! Keep m fixed but increase order:
msq = m*m
do np = m+2, n
  dn = dble(2*np)
  fac = (dn+1.d0)/((dn-3.d0)*dble(np*np-msq))
  nm1 = (np-1)**2
  alp = sqrt(fac*dble(4*nm1-1))
  bet = sqrt(fac*dble(nm1-msq))
  p = alp*z*pmin1-bet*pmin2
  pmin2 = Pmin1
  pmin1 = p
enddo

return
end subroutine generate_spherical

end module spherical
